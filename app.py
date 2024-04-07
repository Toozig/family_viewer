from shiny.express import input, render, ui
from shiny import reactive, req
from family_viewer import  currentState
from viewer_function import get_peak_plot, MAX_SCORE_TFBS, SHOW_SEQ


FAMILY = 'family'
PEAK = 'peak'
SOURCE = 'source'
CHECKBOX = 'checkbox'
N_LINES = 'n_lines'
SCORE_THRESHOLD = 'score_threshold'


CHECKBOX_OPT = [MAX_SCORE_TFBS,
                    SHOW_SEQ]


app_stats = currentState()

family_list = currentState.get_family_list()
source_list = currentState.get_source_list()

# todo: fix this - nreed to init the object and then get the list
peak_list = app_stats.get_peak_list()
min_val, max_val = app_stats.get_threshold_min_max()





ui.page_opts(title="Family Viewer", fillable=True)
with ui.sidebar(width='25vh'):
    # "Sidebar (input)"
    ui.input_selectize(FAMILY, "Family", family_list)
    ui.input_selectize(PEAK, 'Peak',peak_list)
    ui.input_selectize(SOURCE, 'Source',source_list)
    ui.input_slider(SCORE_THRESHOLD, 'Score Threshold',min=min_val, max=max_val, step=1, value=400)
    ui.input_checkbox_group(CHECKBOX, 'Options', CHECKBOX_OPT)
    #ui.input_numeric(N_LINES, 'Number of lines', 2)


# this is how to change the peak list based on the family selection
@reactive.effect
def change_peak_list():
    app_stats.set_family_id(input.family())
    choices = app_stats.get_peak_list()
    ui.update_selectize(PEAK, choices=choices)

# this is to update the slider based on the source and peak selection
@reactive.effect
def change_score_threshold():
    print(f'changing score threshold for {input.source()} and {input.peak()}')
    min_score, max_score = app_stats.get_threshold_min_max()
    min_score, max_score = int(min_score), int(max_score)
    ui.update_slider(SCORE_THRESHOLD, min=min_score, max=max_score, value=min(400, max_score))

@reactive.calc
def update_family_details():
    print(f'updating family details for {input.family()}')
    cur_fam_df =  app_stats.get_family_metadata()
    return cur_fam_df

@reactive.calc
def update_peak_data():
    print(f'updating peak data for {input.peak()}')
    app_stats.set_peak_id(input.peak())
    cur_peak_df =  app_stats.get_peak_data()
    # cur_peak_df =  get_peak_data(input.peak())
    return cur_peak_df

@reactive.calc
def update_variant_list():
    print(f'updating variant list for {input.peak()}, {input.family()}')
    cur_var_df =  app_stats.get_variant_df()
    return cur_var_df

@reactive.effect
def set_source():
    app_stats.set_source(input.source())


def get_peak_plot(peak, source, score_threshold, family, checkbox):
        print(f'getting plot')
        print(f'peak: {peak}, source: {source}, score_threshold: {score_threshold}, family: {family}')
        app_stats.get_peak_plot(checkbox)

with ui.layout_columns(col_widths=[8,4, 12],height='20vh'):
    with ui.card(full_screen=True):

        ui.card_header("family data")

        @render.data_frame
        def family_details():
            return render.DataGrid(update_family_details(),  row_selection_mode="multiple")
        

    with ui.card(full_screen=True):
        ui.card_header("peak data")
        @render.data_frame
        def peak_details():
            df = update_peak_data()[['INTERVAL_ID','CHROM','from','to','length']]
            return render.DataGrid(df,  row_selection_mode="single")
        @render.data_frame
        def peak_details2():
            df = update_peak_data()[['DSDgenes_1_5mb','n_probands','n_healthy']]
            return render.DataGrid(df,  row_selection_mode="single")

with ui.card(full_screen=True):
    # ui.card_header("peak data")
  
    @render.plot()
    def track_plot():
        print('updating plot')
        get_peak_plot(input.peak(), input.source(), input.score_threshold(), input.family(), input.checkbox())

with ui.layout_columns(col_widths=[8, 4, 12],height='20vh'):
    with ui.card():
        ui.card_header("Peak variants")
        @render.data_frame
        def variant_list():
            return render.DataGrid(update_variant_list(),  row_selection_mode="single")


    @reactive.calc
    def update_variant_information():
        # req(input.variant_list_selected_rows())
        print(f"change description for{input.peak(),input.family()}")
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
