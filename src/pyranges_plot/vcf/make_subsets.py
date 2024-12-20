import plotly.graph_objects as go

def make_scatter(p, #vcf_df, 
                 x: str = 'Start',
                 y: str | None=None, # --> count
                 color_by: str | None=None, 
                 size_by: str | None=None,
                 title: str | None=None,
                 title_size: int | None=None,
                 title_color: str | None=None,
                 y_axis_len: int | None=None,
                 y_space: int | None=None,
                 ):
    """
    Create a Scatter plot from a VCF-like DataFrame using Plotly.

    This function generates a scatter plot for visualizing genomic variants or other 
    data points based on the provided DataFrame. It allows customization of axes, 
    marker sizes, colors, and plot titles.

    Parameters
    ----------
        p (pd.DataFrame): 
            Input DataFrame containing the genomic data with columns for 
            start and end positions (e.g., 'Start' and 'End').
        x (str, optional): 
            The column name to use for the x-axis. Defaults to 'Start'.
        y (str): 
            The column name to use for the y-axis.
        color_by (str | None, optional): 
            The column name to use for coloring the markers. If specified, it 
            aggregates unique positions based on this column. Defaults to None.
        size_by (str | None, optional): 
            The column name to use for setting the marker sizes. If specified, it 
            aggregates unique positions based on this column. Defaults to None.
        title (str | None, optional): 
            The title of the plot. Defaults to None.
        title_size (int | None, optional): 
            The font size of the plot title. Applicable only if `title` is specified.
            Defaults to None.
        title_color (str | None, optional): 
            The color of the plot title. Applicable only if `title` is specified.
            Defaults to None.

    Returns
    -------
        Union[go.Scatter, tuple]: 
            - Returns a tuple with the `go.Scatter` object 
              and a dictionary containing title customization options

    Raises
    ------
        ValueError: 
            If `x`, `y`, `color_by`, or `size_by` columns are not found in the 
            input DataFrame.

    Examples
    --------

    """
    # Validate y in input
    if not y:
        raise ValueError(f"The parameter 'y' is required and must be a str to run this function.")
    
    # Validate the x column
    if x not in p.columns:
        raise ValueError(f"The column '{x}' does not exist in the DataFrame.")

    # Handle y-axis logic
    if y not in p.columns:
        raise ValueError(f"The column '{y}' does not exist in the DataFrame.")

     # Add coloring logic if `color_by` is provided
    if color_by:
        if color_by not in p.columns:
            raise ValueError(f"The column '{color_by}' does not exist in the DataFrame.")
        
        # Aggregate color information for unique positions
        color_values = p.groupby(['Start', 'End'])[color_by].first().reset_index()[color_by]
        color_values = color_values.astype('category').cat.codes  # Convert to numeric if categorical
    else:
        color_values = 'blue'  # Default color for all points

    # Handle sizing logic
    if size_by:
        if size_by not in p.columns:
            raise ValueError(f"The column '{size_by}' does not exist in the DataFrame.")
        
        # Aggregate size information for unique positions
        size_values = p.groupby(['Start', 'End'])[size_by].first().reset_index()[size_by]
        size_values = size_values.astype(float)  # Ensure the size column is numeric
    else:
        size_values = 8  # Default marker size

    # Create a scatter plot
    scatter = go.Scatter(
        x=p[x],  # X-axis: Start positions                 #### x 
        y=p[y],  # Y-axis: Counts of transcripts           #### y or __count__
        mode='markers',  # Display points as markers
        marker=dict(
            size=size_values,
            color=color_values,  # Assign color values
            colorscale='Viridis' if color_by else None,  # Use a colormap if coloring
        ),
        hovertemplate='<b>Position:</b> %{x}<br><b>Count:</b> %{y}<extra></extra>'   
    )

    
    custom = {
        'title': title if title else y
    }
    # Defining optional parameters for customisation
    optional_params = {
    'title_size': title_size,
    'title_color': title_color,
    'y_axis_len': y_axis_len,
    'y_space': y_space,
    }

    custom.update({key: value for key, value in optional_params.items() if value is not None})

    return (scatter,custom)