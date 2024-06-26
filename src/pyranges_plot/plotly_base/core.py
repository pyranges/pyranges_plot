from dash import Dash, dcc, html, Input, Output
import dash_bootstrap_components as dbc


def coord2percent(fig, trace, X0, X1):
    """Provides the plot percentage length from the points given. Plotly friendly"""

    if trace == 1:
        trace = ""
    else:
        trace = str(trace)

    x_min, x_max = fig["layout"]["xaxis" + trace]["range"]
    x_rang = x_max - x_min

    percent_size = float(((X1 - X0) / x_rang))

    return percent_size


def percent2coord(fig, trace, x_percent):
    """Provides the coordinates distance from the plot percentage given. Plotly friendly"""

    if trace == 1:
        trace = ""
    else:
        trace = str(trace)

    x_min, x_max = fig["layout"]["xaxis" + trace]["range"]
    x_rang = x_max - x_min

    percent_coord = float(x_percent * x_rang)

    return percent_coord


# Plotly - Function to initialize Dash app layout and callbacks
def initialize_dash_app(fig, max_shown):
    app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

    # create alert and graph components
    subdf_alert = dbc.Alert(
        "The provided data contains more genes than the ones plotted.",
        id="alert-subset",
        color="warning",
        dismissable=True,
        is_open=False,
    )

    uncol_alert = dbc.Alert(
        "Some genes do not have a color assigned so they are colored in black.",
        id="alert-uncolored",
        color="warning",
        dismissable=True,
        is_open=False,
    )

    iter_alert = dbc.Alert(
        "The genes are colored by iterating over the given color list.",
        id="alert-iteration",
        color="warning",
        dismissable=True,
        is_open=False,
    )

    gr = dcc.Graph(id="genes-plot", figure=fig, style={"height": "800px"})

    # define layout
    app.layout = html.Div(
        [dbc.Row([subdf_alert, uncol_alert, iter_alert, gr], justify="around")]
    )

    # callback function
    @app.callback(
        Output("alert-subset", "is_open"),
        Input("genes-plot", "figure"),
    )
    def show_subs_warning(grfig):
        if grfig["data"][0]["customdata"][0] == "no warnings":
            return False
        warn = int(grfig["data"][0]["customdata"][0])
        if warn:
            return True

    @app.callback(
        Output("alert-uncolored", "is_open"),
        Input("genes-plot", "figure"),
    )
    def show_uncol_warning(grfig):
        if grfig["data"][0]["customdata"][0] == "no warnings":
            return False
        sign = int(grfig["data"][0]["customdata"][1])
        if sign == 91124:
            return True

    @app.callback(
        Output("alert-iteration", "is_open"),
        Input("genes-plot", "figure"),
    )
    def show_iter_warning(grfig):
        if grfig["data"][0]["customdata"][0] == "no warnings":
            return False
        sign = int(grfig["data"][0]["customdata"][2])
        if sign == 91321:
            return True

    return app
