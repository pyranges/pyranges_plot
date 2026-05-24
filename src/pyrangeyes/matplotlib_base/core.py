import os
import ipywidgets as widgets
from IPython.display import display

if os.environ.get("DISPLAY"):
    try:
        import tkinter as tk
    except ImportError:
        tk = None
else:
    tk = None


def coord2percent(ax, X0, X1):
    """Provides the plot percentage length from the points given. Matplotlib friendly"""

    x_min, x_max = ax.get_xlim()
    x_rang = x_max - x_min
    percent_size = float(((X1 - X0) / x_rang))

    return percent_size


def percent2coord(ax, x_percent):
    """Provides the coordinates distance from the plot percentage given. Matplotlib friendly"""

    x_min, x_max = ax.get_xlim()
    x_rang = x_max - x_min
    percent_coord = float(x_percent * x_rang)

    return percent_coord


def running_in_jupyter():
    try:
        shell = get_ipython().__class__.__name__
        return shell == "ZMQInteractiveShell"
    except NameError:
        return False


if running_in_jupyter():

    def plt_popup_warning(txt, bkg="#1f1f1f", txtcol="white", botcol="#D6AA00"):
        # Widget Label i Botó

        label = widgets.HTML(
            value=f'<div style="color:{txtcol}; background-color:{bkg}; padding:10px; font-family:Sans; font-size:15px;">{txt}</div>'
        )
        button = widgets.Button(
            description="Got it",
            style={"button_color": botcol, "font_weight": "bold"},
            layout=widgets.Layout(margin="10px 0 0 0"),
        )

        box = widgets.HBox([label, button])
        display(box)

        def _on_click(b):
            try:
                box.close()
            except Exception:
                pass

        button.on_click(_on_click)

else:

    def plt_popup_warning(txt, bkg="#1f1f1f", txtcol="white", botcol="#D6AA00"):
        warn = tk.Tk()
        warn.wm_title("Warning!")
        warn.configure(background=bkg)

        label = tk.Label(warn, label=txt, font=("Sans", 15), fg=txtcol, bg=bkg)
        label.pack(side="top", anchor="center", pady=10)

        bot = tk.Button(
            warn, label="Got it", command=warn.destroy, fg="black", bg=botcol
        )
        bot.pack(pady=10)

        warn.wait_window()


def make_annotation(item, fig, ax, geneinfo, tag_background):
    """Create annotation for a given plot item."""

    # create annotation and make it not visible
    annotation = ax.annotate(
        "",
        xy=(0, 0),
        xytext=(20, 20),
        textcoords="offset points",
        bbox=dict(
            boxstyle="round",
            edgecolor=tag_background,
            facecolor=tag_background,
        ),
        arrowprops=dict(arrowstyle="->"),
        color="white",
    )
    annotation.set_visible(False)

    # make annotation visible when over the gene line
    def on_hover(event):
        visible = annotation.get_visible()
        contains_item, _ = item.contains(event)  # Check if mouse is over the gene line
        if contains_item:
            annotation.set_text(geneinfo)
            annotation.xy = (event.xdata, event.ydata)
            annotation.set_visible(True)
            fig.canvas.draw()
        elif visible:
            annotation.set_visible(False)
            fig.canvas.draw()

    fig.canvas.mpl_connect("motion_notify_event", on_hover)


def rgb_string_to_tuple(rgb_string):
    # Store the numbers
    rgb_string = rgb_string.lstrip("rgb(").rstrip(")")

    # Split the string by commas and convert to floats
    r, g, b = map(int, rgb_string.split(","))

    # Convert to floats in the range [0, 1]
    return (r / 255.0, g / 255.0, b / 255.0)
