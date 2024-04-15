from shiny.express import input, render, ui
from shiny import reactive, req
from family_viewer import  currentState
from viewer_function import create_legend

FAMILY = 'family'
PEAK = 'peak_list'



app_stats = currentState()

family_list = currentState.get_family_list()

peak_list = app_stats.get_peak_list()
min_val, max_val = app_stats.get_threshold_min_max()





ui.page_opts(title="Family Viewer", full_width=True)


update_plot = reactive.value(False)
var_df = reactive.value(app_stats.var_df)
peak = reactive.value(app_stats.var_df.to_dict('records')[0])


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


ui.tags.style('.small_font { font-size: 1.3vh}')
with ui.layout_columns(col_widths=[1,5,6,13],height='15vh', class_='small_font'):

    
    ui.input_select(FAMILY, "Family", family_list)

    @render.data_frame
    def family_details():
        return render.DataGrid(update_family_details(),  row_selection_mode="multiple")
    
    with ui.card():
        @render.data_frame
        def peak_details():
            df = update_peak_data()[['INTERVAL_ID','CHROM','from','to','length']]
            return render.DataGrid(df,  row_selection_mode="single")
        @render.data_frame
        def peak_details2():
                df = update_peak_data()[['distance_from_nearest_DSD_TSS','n_probands','n_healthy']]
                return render.DataGrid(df,  row_selection_mode="single")
    
with ui.card(class_='small_font'):
    ui.card_header("family variants")

    @render.data_frame()
    def family_variants():
        print(f'updating family variants for {input.family()}')
        app_stats.set_family_id(input.family())
        df = app_stats.get_view_all_variants()
        df = df.reset_index()
        df = df.rename(columns={'index':'variant_id'})
        df.REF = df.REF.apply(lambda x: x[:7])
        df.ALT = df.ALT.apply(lambda x: x[:7])
        var_df.set(app_stats.var_df.iloc[df.index,:])
        return render.DataGrid(df,  row_selection_mode="single", height='20vh', width='700vh')


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



with ui.card(full_screen=True):
        
    with ui.layout_column_wrap(width= 1 /7, gap='1vh',
                               fixed_width='TRUE',
                                height='5vh',
                                class_="d-flex justify-content-between align-items-center"):
        ui.input_action_button("upstream", "<<<<<", width='100%', height='100%' )
        ui.input_numeric('add_bp', '    bp', 500)
        ui.input_action_button("downstream", ">>>>>", width='100%', height='100%')
        with ui.layout_column_wrap(width= 1 /2, height='100%', gap='4px'):
            ui.input_switch("zoom", "zoom", True) 
            @render.ui
            def zoom():
                zoom_status = "Out" if input.zoom() else "In"
                return zoom_status
        
        ui.input_action_button("reset", "reset view")

        @render.plot()
        def legend_plot():
            return create_legend()
       

    @render.plot()  
    def track_plot():  
        print('updating plot')
        update_plot.get()
        print(f'{input.downstream()}')
        print(f'{input.upstream()}')
        input.reset()
        get_track_plot(peak.get())




@render.ui
def var_id_text():
    return f'all variants in {peak.get()["INTERVAL_ID"]}'

with ui.layout_columns(col_widths=[8,4, 12],height='25vh',class_='small_font'):

    # with ui.card():
        # ui.card_header("peak data")
        # @render.data_frame
        # def peak_details():
        #     df = update_peak_data()[['INTERVAL_ID','CHROM','from','to','length']]
        #     return render.DataGrid(df,  row_selection_mode="single")
        
        # @render.data_frame
        # def peak_details2():
        #     df = update_peak_data()[['distance_from_nearest_DSD_TSS','n_probands','n_healthy']]
        #     return render.DataGrid(df,  row_selection_mode="single")

    with ui.layout_column_wrap(width= "8/1",heights_equal='row', fill=False):
        ui.input_checkbox('focus','focus on click', False)

        @render.data_frame()
        def all_peak_variants():
            print(f'updating peaks variants for {peak.get()["INTERVAL_ID"]}')
            df = app_stats.get_all_peaks_df()

            return render.DataGrid(df,  row_selection_mode="single", height='25vh', width='700vh')

    with ui.layout_column_wrap(width= "4/1",):
        @render.ui
        def var_information():
            peak.get()
            df = app_stats.get_all_peak_vars( to_filter = False)
            index_list = input.all_peak_variants_selected_rows()
            index = 0 if len(index_list) == 0 else index_list[0]
            string = app_stats.get_variant_info(index,df)
            print(string)   
            return ui.markdown(string)
        # ui.markdown("** hey! **")

        @reactive.effect
        @reactive.event( input.all_peak_variants_selected_rows)
        def change_focus():
            req(input.all_peak_variants_selected_rows())
            if input.focus():
                index = input.all_peak_variants_selected_rows()[0]
                print(f'changing focus to {index}')
                variant = app_stats.get_all_peak_vars( to_filter = False).to_dict('records')[index]
                variant['from'] = variant['POS'] - 250
                variant['to'] = variant['POS'] + 250
                peak.set(variant)


    # with ui.card():
    #     ui.card_header("Peak variants")
    #     @render.data_frame
    #     def variant_list():
    #         print(f'updating variant list for {peak.get()["INTERVAL_ID"]}')
    #         df = app_stats.get_variant_df_by_id(peak.get()['INTERVAL_ID'])
    #         # df = df[df.INTERVAL_ID == peak.get()['INTERVAL_ID']]
            
    #         df = df.reset_index()
    #         df['index'] = df['index'] + 1
    #         return render.DataGrid(df,  row_selection_mode="single")


    # @reactive.calc
    # def update_variant_information(index, df=None):
    #     # req(input.variant_list_selected_rows())
    #     print(f"change description for{peak.get()['INTERVAL_ID'],input.family()}")
    #     string = app_stats.get_variant_info(index,df)
    #     return string
