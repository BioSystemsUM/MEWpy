import matplotlib
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class Plot:

    def __init__(self,
                 plot_title='Pareto Aproximation',
                 reference_front=None,
                 reference_point=None,
                 axis_labels=None):
        """
        :param plot_title: Title of the graph.
        :param axis_labels: List of axis labels.
        :param reference_point: Reference point (e.g., [0.4, 1.2]).
        :param reference_front: Reference Pareto front (if any) as solutions.
        """
        self.plot_title = plot_title
        self.axis_labels = axis_labels
        self.reference_point = reference_point
        self.reference_front = reference_front
        self.dimension = None

    @staticmethod
    def get_points(solutions):
        """ Get points for each solution of the front.

        :param solutions: List of solutions where each solution is a list of fitness values
        :return: Pandas dataframe with one column for each objective and one row for each solution.
        """
        p = [s.fitness for s in solutions]
        points = pd.DataFrame(p)
        return points, points.shape[1]

    def plot(self, front, label='', normalize=False, filename=None, format='eps'):
        """ Plot any arbitrary number of fronts in 2D, 3D or p-coords.

        :param front: Pareto front or a list of them.
        :param label: Pareto front title or a list of them.
        :param normalize: If True, normalize data (for p-coords).
        :param filename: Output filename.
        :param format: Output file format.
        """
        if not isinstance(front[0][0], float):
            front = [front]

        if not isinstance(label, list):
            label = [label]

        if len(front) != len(label):
            raise Exception('Number of fronts and labels must be the same')

        dimension = len(front[0][0])

        if dimension == 2:
            self.two_dim(front, label, filename, format)
        elif dimension == 3:
            self.three_dim(front, label, filename, format)
        else:
            self.pcoords(front, normalize, filename, format)

    def two_dim(self, fronts, labels=None, filename=None, format='eps'):
        """ Plot any arbitrary number of fronts in 2D.

        :param fronts: List of fronts (containing solutions).
        :param labels: List of fronts title (if any).
        :param filename: Output filename.
        """
        n = int(np.ceil(np.sqrt(len(fronts))))
        fig = plt.figure()
        fig.suptitle(self.plot_title, fontsize=16)

        reference = None
        if self.reference_front:
            reference, _ = self.get_points(self.reference_front)

        for i, _ in enumerate(fronts):
            points, _ = self.get_points(fronts[i])

            ax = fig.add_subplot(n, n, i + 1)
            points.plot(kind='scatter', x=0, y=1, ax=ax, s=10, color='#236FA4', alpha=1.0)

            if labels:
                ax.set_title(labels[i])

            if self.reference_front:
                reference.plot(x=0, y=1, ax=ax, color='k', legend=False)

            if self.reference_point:
                for point in self.reference_point:
                    plt.plot([point[0]], [point[1]], marker='o', markersize=5, color='r')
                    plt.axvline(x=point[0], color='r', linestyle=':')
                    plt.axhline(y=point[1], color='r', linestyle=':')

            if self.axis_labels:
                plt.xlabel(self.axis_labels[0])
                plt.ylabel(self.axis_labels[1])

        if filename:
            plt.savefig(filename + '.' + format, format=format, dpi=200)

        plt.show()
        plt.close(fig)

    def three_dim(self, fronts, labels=None, filename=None, format='eps'):
        """ Plot any arbitrary number of fronts in 3D.

        :param fronts: List of fronts (containing solutions).
        :param labels: List of fronts title (if any).
        :param filename: Output filename.
        """
        n = int(np.ceil(np.sqrt(len(fronts))))
        fig = plt.figure()
        fig.suptitle(self.plot_title, fontsize=16)

        for i, _ in enumerate(fronts):
            ax = fig.add_subplot(n, n, i + 1, projection='3d')
            ax.scatter([s.objectives[0] for s in fronts[i]],
                       [s.objectives[1] for s in fronts[i]],
                       [s.objectives[2] for s in fronts[i]])

            if labels:
                ax.set_title(labels[i])

            if self.reference_front:
                ax.scatter([s.objectives[0] for s in self.reference_front],
                           [s.objectives[1] for s in self.reference_front],
                           [s.objectives[2] for s in self.reference_front])

            if self.reference_point:
                # todo
                pass

            ax.relim()
            ax.autoscale_view(True, True, True)
            ax.view_init(elev=30.0, azim=15.0)
            ax.locator_params(nbins=4)

        if filename:
            plt.savefig(filename + '.' + format, format=format, dpi=1000)

        plt.show()
        plt.close(fig)

    def pcoords(self, fronts, normalize=False, filename=None, format='eps'):
        """ Plot any arbitrary number of fronts in parallel coordinates.

        :param fronts: List of fronts (containing solutions).
        :param filename: Output filename.
        """
        n = int(np.ceil(np.sqrt(len(fronts))))
        fig = plt.figure()
        fig.suptitle(self.plot_title, fontsize=16)

        for i, _ in enumerate(fronts):
            points, _ = self.get_points(fronts[i])

            if normalize:
                points = (points - points.min()) / (points.max() - points.min())

            ax = fig.add_subplot(n, n, i + 1)
            pd.plotting.parallel_coordinates(points, 0, ax=ax)

            ax.get_legend().remove()

            if self.axis_labels:
                ax.set_xticklabels(self.axis_labels)

        if filename:
            plt.savefig(filename + '.' + format, format=format, dpi=1000)

        plt.show()
        plt.close(fig)


class StreamingPlot:

    def __init__(self,
                 plot_title='Pareto Approximation',
                 reference_front=None,
                 reference_point=None,
                 axis_labels=None):
        """
        :param plot_title: Title of the graph.
        :param axis_labels: List of axis labels.
        :param reference_point: Reference point (e.g., [0.4, 1.2]).
        :param reference_front: Reference Pareto front (if any) as solutions.
        """
        self.plot_title = plot_title
        self.axis_labels = axis_labels

        if reference_point and not isinstance(reference_point[0], list):
            reference_point = [reference_point]

        self.reference_point = reference_point
        self.reference_front = reference_front
        self.dimension = None

        import warnings
        warnings.filterwarnings("ignore", ".*GUI is implemented.*")

        self.fig, self.ax = plt.subplots()
        self.sc = None
        self.scf = None
        self.axis = None

    def plot(self, front, dominated=None):
        # Get data
        points, dimension = Plot.get_points(front)

        # Create an empty figure
        self.create_layout(dimension)

        # If any reference point, plot
        if self.reference_point:
            for point in self.reference_point:
                self.scp, = self.ax.plot(*[[p] for p in point], c='r', ls='None', marker='*', markersize=3)

        # Plot data
        self.sc, = self.ax.plot(*[points[column].tolist() for column in points.columns.values],
                                ls='None', marker='o', markersize=4)

        # dominated points
        # if dominated:
        #    rpoints, _ = Plot.get_points(dominated)
        #    self.scf, = self.ax.plot(*[[p] for p in rpoints],
        #                             c='g', ls='None', marker='o', markersize=1)

        # Show plot
        plt.show(block=False)

    def update(self, front, dominated=None, reference_point=None, text=None):
        if self.sc is None:
            raise Exception('Figure is none')

        points, dimension = Plot.get_points(front)

        # if dominated:
        #    dpoints = Plot.get_points(dominated)
        #    self.scf.set_data(dpoints[0], dpoints[1])
        #    if dimension == 3:
        #        self.scf.set_3d_properties(dpoints[2])

        # Replace with new points
        self.sc.set_data(points[0], points[1])

        if dimension == 3:
            self.sc.set_3d_properties(points[2])

        # If any new reference point, plot
        if reference_point:
            self.scp.set_data([p[0] for p in reference_point], [p[1] for p in reference_point])

        #if text is not None:
        #    self.fig.title(text, fontsize=10)

        # Re-align the axis
        self.ax.relim()
        self.ax.autoscale_view(True, True, True)

        try:
            # self.fig.canvas.draw()
            self.fig.canvas.flush_events()
        except KeyboardInterrupt:
            pass

        pause(0.01)

    def create_layout(self, dimension):
        self.fig.canvas.set_window_title(self.plot_title)
        self.fig.suptitle(self.plot_title, fontsize=16)

        if dimension == 2:
            # Stylize axis
            self.ax.spines['top'].set_visible(False)
            self.ax.spines['right'].set_visible(False)
            self.ax.get_xaxis().tick_bottom()
            self.ax.get_yaxis().tick_left()
            if self.axis_labels:
                plt.xlabel(self.axis_labels[0])
                plt.ylabel(self.axis_labels[1])
        elif dimension == 3:
            self.ax = Axes3D(self.fig)
            self.ax.autoscale(enable=True, axis='both')
            if self.axis_labels:
                self.ax.set_xlabel(self.axis_labels[0])
                self.ax.set_ylabel(self.axis_labels[1])
                self.ax.set_zlabel(self.axis_labels[2])
        else:
            raise Exception('Dimension must be either 2 or 3')

        self.ax.set_autoscale_on(True)
        self.ax.autoscale_view(True, True, True)

        # Style options
        self.ax.grid(color='#f0f0f5', linestyle='-', linewidth=0.5, alpha=0.5)


def pause(interval):
    backend = plt.rcParams['backend']

    if backend in matplotlib.rcsetup.interactive_bk:
        figManager = matplotlib._pylab_helpers.Gcf.get_active()
        if figManager is not None:
            canvas = figManager.canvas
            if canvas.figure.stale:
                canvas.draw()
            canvas.start_event_loop(interval)
            return
