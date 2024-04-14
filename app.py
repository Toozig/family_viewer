from shiny.express import input, render, ui
from shiny import reactive, req
from family_viewer import  currentState
from viewer_function import  MAX_SCORE_TFBS, SHOW_SEQ

FAMILY = 'family'
PEAK = 'peak_list'



app_stats = currentState()

family_list = currentState.get_family_list()

peak_list = app_stats.get_peak_list()
min_val, max_val = app_stats.get_threshold_min_max()





ui.page_opts(title="Family Viewer", fillable=True)


update_plot = reactive.value(False)
peak = reactive.value(app_stats.peak_id)
var_df = reactive.value(app_stats.var_df)

@reactive.effect
@reactive.event(input.downstream)
def add_downstream():
    print(input.downstream())
    add_value = input.add_bp() *  (1 if input.zoom() else -1)
    app_stats.add_downstream(add_value)
    print(f'adding downstream {add_value}')

@reactive.effect
@reactive.event(input.upstream)
def add_upstream():
    print(input.upstream())
    add_value = input.add_bp() * ( 1 if input.zoom() else -1)
    app_stats.add_upstream(add_value)


@reactive.effect
@reactive.event(input.reset)
def reset_view():
    print(input.reset())
    app_stats.reset_view()

@reactive.calc
def update_family_details():
    print(f'updating family details for {input.family()}')
    cur_fam_df =  app_stats.get_family_metadata()
    return cur_fam_df

@reactive.calc
def update_peak_data():
    # print(f'updating peak data for {peak.get()["INTERVAL_ID"]}')
    app_stats.set_peak_id(peak.get()['INTERVAL_ID'])
    cur_peak_df =  app_stats.get_peak_data()
    app_stats.reset_view()

    return cur_peak_df

@reactive.calc
def update_variant_list():
    # print(f'updating variant list for {peak.get()['INTERVAL_ID']}, {input.family()}')
    cur_var_df =  app_stats.get_variant_df()
    return cur_var_df


def get_track_plot(peak):
        print(f'getting plot')
        print(f'peak: {peak}, source: family:')
        app_stats.get_track_plot(peak)


with ui.layout_columns(col_widths=[8,4, 12],height='20vh'):

    with ui.card():
        with ui.layout_columns():
            ui.input_select(FAMILY, "Family", family_list)
            # ui.input_select(PEAK, 'Peak',peak_list)
        ui.card_header("family variants")

        
        @render.data_frame
        def family_variants():
            print(f'updating family variants for {input.family()}')
            app_stats.set_family_id(input.family())
            df = app_stats.get_view_all_variants()
            df = df.reset_index()
            df = df.rename(columns={'index':'variant_id'})
            var_df.set(app_stats.var_df.iloc[df.index,:])
            return render.DataGrid(df,  row_selection_mode="single")


        @reactive.effect
        @reactive.event(input.family_variants_selected_rows)
        def change_track_by_family_var():
            selected_row = input.family_variants_selected_rows()
            print(f'changing track by family variant {selected_row}')
            index = 0 if len(selected_row) == 0 else selected_row[0]
            print(f'changing focus to {index}')
            variant = var_df.get().to_dict('records')[index]
            # print(variant)
            peak.set(variant)
            print(f'changing focus to {variant["INTERVAL_ID"]}')

    with ui.card():
        ui.card_header("family data")
        @render.data_frame
        def family_details():
            return render.DataGrid(update_family_details(),  row_selection_mode="multiple")
        


with ui.layout_columns(col_widths=[1,11, 12],height='25vh'):

    with ui.card(height='10vh'):
        ui.input_numeric('add_bp', 'see more bp', 500)
        with ui.layout_columns():
            ui.input_action_button("upstream", "upstream" )
            ui.input_action_button("downstream", "downstream")
        ui.input_switch("zoom", "zoom", True) 

        @render.ui
        def value():
            zoom_status = "Out" if input.zoom() else "In"
            return zoom_status
            
        ui.input_action_button("reset", "reset view")
        ui.input_checkbox('focus','focus on click', True)

    with ui.card(full_screen=True):
        # ui.card_header("peak data")
        @render.plot()
        def track_plot():
            print('updating plot')
            update_plot.get()
            print(f'{input.downstream()}')
            print(f'{input.upstream()}')
            input.reset()
            get_track_plot(peak.get())

with ui.layout_columns(col_widths=[4,4, 4, 12],height='15vh'):

    with ui.card():
        ui.card_header("peak data")
        @render.data_frame
        def peak_details():
            df = update_peak_data()[['INTERVAL_ID','CHROM','from','to','length']]
            return render.DataGrid(df,  row_selection_mode="single")
        
        @render.data_frame
        def peak_details2():
            df = update_peak_data()[['distance_from_nearest_DSD_TSS','n_probands','n_healthy']]
            return render.DataGrid(df,  row_selection_mode="single")
        

    with ui.card():
        ui.card_header("Peak variants")
        @render.data_frame
        def variant_list():
            print(f'updating variant list for {peak.get()["INTERVAL_ID"]}')
            df = app_stats.get_variant_df_by_id(peak.get()['INTERVAL_ID'])
            # df = df[df.INTERVAL_ID == peak.get()['INTERVAL_ID']]
            
            df = df.reset_index()
            df['index'] = df['index'] + 1
            return render.DataGrid(df,  row_selection_mode="single")


    @reactive.calc
    def update_variant_information():
        # req(input.variant_list_selected_rows())
        print(f"change description for{peak.get()['INTERVAL_ID'],input.family()}")
        index = 0 if len(input.variant_list_selected_rows()) == 0 else input.variant_list_selected_rows()[0]
        string = app_stats.get_variant_info(index)
        return string

    with ui.card(width='15vh'):

        ui.card_header("variant information")
        @render.ui
        def var_information():
            string = update_variant_information()
            print(string)   
            return ui.markdown(string)
        # ui.markdown("** hey! **")

    @reactive.effect
    @reactive.event(input.variant_list_selected_rows)
    def change_focus():
        req(input.variant_list_selected_rows())
        if input.focus():
            index = input.variant_list_selected_rows()[0]
            print(f'changing focus to {index}')
            variant = update_variant_list().copy().to_dict('records')[index]
            variant['from'] = variant['POS'] - 250
            variant['to'] = variant['POS'] + 250
            peak.set(variant)

