import json
import argparse
import os
import pandas as pd
from bokeh.io import output_file, show, save
from bokeh.models import HoverTool
from bokeh.plotting import figure, ColumnDataSource
from bokeh.palettes import Category20c

"""If importing this module into a Jupyter Notebook, import the following:
from ipywidgets import interact
from IPython.display import display, HTML
from bokeh.io import push_notebook, output_notebook
output_notebook() # for plotting Bokeh in notebook
"""


class Distribution(object):
    """A structure to construct and hold the fragment length distribution.

    Attributes:
        file_path: The path to the JSON file holding the count information
    """

    def __init__(self, file_path):
        """Returns a new Distribution object"""
        self._path = file_path
        self._d = self.open_json()
        self._name = '_'.join(self._path.split('/')[-1].split('_')[:2])
        self._lengths = self.mapping.keys()
        self._counts = self.mapping.values()
        self._fragments = float(sum(self.counts))
        self._densities = {l: c / self.fragments
                           for l, c in zip(self.lengths, self.counts)}

    @property
    def mapping(self):
        # type: () -> dict[int, float]
        """Returns a dictionary that maps fragment length to it's count"""
        return self._d

    @property
    def name(self):
        # type: () -> str
        """Returns the name of the sample this Distribution is for"""
        return self._name

    @property
    def lengths(self):
        # type: () -> list[int]
        """Returns a list of the fragment lengths"""
        return self._lengths

    @property
    def counts(self):
        # type: () -> list[float]
        """Returns a list of counts. Index is the same as self.lengths"""
        return self._counts

    @property
    def fragments(self):
        # type: () -> float
        """Returns the total number of fragments in sample"""
        return self._fragments

    @property
    def mean(self):
        # type: () -> float
        """Calculates the mean fragment length for this sample

        :return: the mean fragment length
        """
        s = sum([
                l * c
                for l, c in zip(self.lengths, self.counts)
            ])
        return s / self.fragments

    @property
    def densities(self):
        # type: () -> dict[int, float]
        """Returns a dictionary that maps fragment length to it's density"""
        return self._densities

    def to_df(self):
        # type: () -> pd.DataFrame
        """Convert the sample dictionary of fragment length distribution 
        densities into a Pandas DataFrame. Column name will be the sample ID.

        :return: Pandas DataFrame of the fragment length distribution for 
        this sample.
        """
        df = pd.DataFrame.from_dict(self.densities, orient='index')
        df.columns = [self.name]
        return df

    def open_json(self):
        # type: () -> dict[int, float]
        """Load in JSON file as dictionary.

        :return: Dictionary of the contents of the JSON file for this sample
        """
        with open(self._path, 'r') as f:
            d = json.load(f)[0]

        d = {int(k): float(v) for k, v in d.items()}
        return d


def join_dfs(dfs):
    # type: (list[pd.DataFrame]) -> pd.DataFrame
    """Joins multiple Pandas DataFrames into one. 
    NaN values will be changed to 0.

    :param dfs: list of Pandas DataFrames for each sample.
    :return: A Pandas DataFrame of all sample's DataFrames joined together.
    """
    df = pd.concat(dfs, axis=1)
    return df.fillna(value=0)


def make_plot(df, y_label, y_scale, plot_name):
    # type: (pd.DataFrame, str, str, str) -> None
    """Function that creates and saves the figure.

    :param df: Pandas DataFrame of the fragment length distributions for all
    samples.
    :param y_label: Label to use for y-axis.
    :param y_scale: y-axis scaling. i.e 'linear', 'log' etc.
    :param plot_name: name for the HTML file.
    """
    num_lines = len(df.columns)
    legends_list = df.columns
    my_pallete = Category20c[num_lines]
    xs = [df.index.values] * num_lines
    ys = [df[name].values for name in df]

    # create the figure skeleton
    p = figure(title='Fragment Length Distribution',
               x_axis_label='Fragment Length (bp)',
               y_axis_label=y_label,
               y_axis_type=y_scale,
               plot_height=600,
               plot_width=960,
               tools='pan, box_zoom, reset, save, wheel_zoom')

    # add a line to the figure for each sample
    for (col, leg, x, y) in zip(my_pallete, legends_list, xs, ys):
        src = create_data_source(x, y, leg)
        p.line('x', 'y', color=col, legend=leg, source=src)

    # Information to be displayed on the hover tooltips
    tooltips = [
        ('Length', '$index'),
        ('Density', '$y'),
        ('Sample', '@sample')
    ]

    # adjusting figure styling
    p.add_tools(HoverTool(tooltips=tooltips, mode='vline'))
    p.legend.location = 'top_right'
    p.legend.click_policy = 'hide'  # clicking group will hide it
    p.legend.label_text_font = 'roboto'
    p.legend.background_fill_alpha = 0

    # write figure to a html file
    output_file(plot_name)
    # show(p)
    save(p)


def create_data_source(x, y, leg):
    # type: (list[int], list[float], list[str]) -> ColumnDataSource
    """Creates the Bokeh data source object that helps with the tooltip.
    
    :param x: list of the fragment lengths.
    :param y: list of the densities.
    :param leg: sample name that will be used in the legend.
    :return: A Bokeh class ColumnDataSource which is a mapping for the figure
    to use when creating the line.
    """
    return ColumnDataSource({
            'x': x,
            'y': y,
            'sample': [leg] * len(x)
        })


def plot_master(file_paths, threshold, plot_name):
    # type: (list[str], int, str) -> None
    """Function that coordinates creating plots.

    :param file_paths: list of file paths for each JSON file.
    :param threshold: do not plot fragments longer than this threshold.
    :param plot_name: name for the HTML file of the figure.
    """
    # make a dataframe of all sample's fragment length distributions
    df = join_dfs([Distribution(f).to_df() for f in file_paths])

    # drop all information for fragments over given length threshold
    df.drop(df.index[threshold:], inplace=True)

    # create the normal and log plots and writes them to file
    make_plot(df, 'Density', 'linear', plot_name + '.html')
    make_plot(df, 'log10 Density', 'log', plot_name + '_log.html')


def main():
    parser = argparse.ArgumentParser(
        description="""Plot fragment length distribution for cell-free
        DNA samples. Will save plots to HTML file and you can open these in 
        a browser for interactive visualisation."""
    )
    parser.add_argument(
        '-i',
        '--input_dir',
        help="""Directory where JSON files containing the fragment length 
        counts are stored.""",
        required=True,
        type=str
    )
    parser.add_argument(
        '-l',
        '--length_threshold',
        help="""Only plot fragments of a length less than the given threshold.
        Default: 1000""",
        default=1000,
        type=int
    )
    parser.add_argument(
        '-o',
        '--output_name',
        help="""Name of output file name. If plotting log as well, output_name
         will have log appended to it. You may also specify a prefixed path for 
          where the output should be stored. Default: fragment_length""",
        default='fragment_length',
        type=str
    )

    args = parser.parse_args()

    # get list of all json file paths
    file_paths = [os.path.join(args.input_dir, f)
                  for f in os.listdir(args.input_dir)
                  if f.endswith("_count.json")]

    plot_master(file_paths, args.length_threshold, args.output_name)


if __name__ == "__main__":
    main()
