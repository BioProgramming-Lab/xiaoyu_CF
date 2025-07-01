import mpl_scatter_density
import matplotlib.pyplot as plt
import cytoflow
import glob
import os.path
import numpy as np

# suppress warnings
# ZE5 has some problems with the data output, which will cause some warnings.
# but the warnings are not important for the data analysis.
# and the data is just 4 order of magnitude larger than it should be.
import warnings
warnings.filterwarnings("ignore")

class xiaoyu_Expr:

    # load fcs files in the working directory
    # usage: xiaoyu_Expr("working_directory", ["merge_directory1", "merge_directory2", ...])
    def __init__(self, working_dir : str, merge_dirs = []):
        self.working_dir = ""
        self.filelist = []
        self.name_locus_mapping = {}
        self.expr = None
        if working_dir != "" and working_dir[-1] != "/":
            working_dir += "/"
        self.working_dir = working_dir
        filelist = glob.glob(working_dir + "*.fcs")
        filelist = sorted(filelist)
        self.filelist = filelist
        print (f"{len(filelist)} fsc files loaded")

        # load other fcs files in the merge_dirs, and merge them to the current experiment by well index.
        # this is useful when you have some additive data to the current experiment.
        if len(merge_dirs) != 0:
            for merge_dir in merge_dirs:
                if merge_dir[-1] != "/":
                    merge_dir += "/"
                filelist_merge = glob.glob(merge_dir + "*.fcs")
                filelist_merge = sorted(filelist_merge)
                self.filelist += filelist_merge
                print (f"other {len(filelist_merge)} fsc files loaded from {merge_dir}. will be merged to the current experiment by well index.")

        # classify samples with name (A1, A2, A3, ..., B1, B2, B3, ...)
        file_names = []
        for i in filelist:
            file_names.append(os.path.splitext(os.path.split(i)[-1])[0])
            file_names[-1] = file_names[-1].split(" ")[0]   # filter out the well name. e.g. "A1 EGFR.fcs" -> "A1"
        for i in range(1, len(file_names)+1):
            self.name_locus_mapping[i] = [file_names[i-1][0], int(file_names[i-1][1:])-1]
        #print (f"mapping of sample name to char and num: {name_locus_mapping}")

        # for merge_dir in merge_dirs: ------------------------------------------------------------- function addition: merge some additive data to the current experiment
            
        # generate experiment instance
        tubes = []
        count = 1
        # for merged data, the sample index (num, char, well) should be the same with the original data, but the sample index (sample) should be different for clarification.
        for i in filelist:
            tube = cytoflow.Tube(file = i, conditions = {"sample" : count, "char" : self.name_locus_mapping[count][0], "num" : str(self.name_locus_mapping[count][1]+1), "well" : self.name_locus_mapping[count][0]+str(self.name_locus_mapping[count][1]+1)})
            tubes.append(tube)
            count += 1
        experiment = cytoflow.ImportOp(conditions = {"sample" : "int", "char" : "category", "num" : "category", "well" : "str"}, tubes = tubes).apply()
        print (f"properties of the experiment: {experiment.data.columns}")
        self.expr = experiment

    # add concetration condition to the experiment
    # usage: add_conc_condition(ori_conc, dilu_dir = "down", dilu_factor = 5)
    # this function will add a new condition to the experiment with the name "conc" and the dtype "float".
    # the concentration of the samples will be assigned according to the dilution direction and dilution factor.
    # the concentration of the first sample will be the ori_conc, and the concentration of the other samples will be calculated by the dilution factor.
    # this sample only works for 96-well plate.
    def add_conc_condition(self, ori_conc, dilu_dir = "down", dilu_factor = 5):
        conc = []
        condition_list = []
        if dilu_dir == "left" or dilu_dir == "right":
            for i in range(12):
                if dilu_dir == "left":
                    conc.append(ori_conc / (dilu_factor ** (12-i)))
                else:
                    conc.append(ori_conc / (dilu_factor ** i))
            for i in self.expr.data["num"]:
                condition_list.append( round(conc[int(i)-1], 18) )
            condition_list = np.array(condition_list)
            self.expr.add_condition(name = "conc", dtype = "float", data = condition_list)
        if dilu_dir == "down" or dilu_dir == "up":
            for i in range(8):
                if dilu_dir == "up":
                    conc.append(ori_conc / (dilu_factor ** (8-i)))
                else:
                    conc.append(ori_conc / (dilu_factor ** i))
            for i in self.expr.data["char"]:
                condition_list.append( round(conc[ord(i)-ord('A')], 18) )
            condition_list = np.array(condition_list)
            self.expr.add_condition(name = "conc", dtype = "float", data = condition_list)

    # add other conditions
    # usage: add_condition_with_sample_array("condition_name", "dtype", [condition1, condition2, condition3, ...])
    # this function will add a new condition to the experiment with the name "condition_name" and the dtype "dtype".
    # by assigning the list of conditions to their corresponding samples.
    # the length of the conditions should be the same with the number of samples in the experiment.
    # if module tabulate is installed, the mapping of the sample index to the new condition value will be shown.
    def add_condition_with_sample_array(self, condition_name, dtype, conditions):
        condition_list = []
        for i in self.expr.data["sample"]:
            condition_list.append( conditions[i-1] )
        condition_list = np.array(condition_list)
        self.expr.add_condition(name = condition_name, dtype = dtype, data = condition_list)
        print (f"condition {condition_name} added to the experiment")
        print ("mapping of the sample index to the new condition value is shown below:")
        try:
            from tabulate import tabulate
            print (tabulate([[i, f"{self.name_locus_mapping[i][0]}{self.name_locus_mapping[i][1]+1}", conditions[i-1]] for i in range(1, len(self.name_locus_mapping)+1)], headers = ["sample", "well", condition_name], tablefmt='fancy_grid'))
        except:
            print ("tabulate module not found. please install tabulate module to show the mapping of the sample index to the new condition value.")

    def sample_index_to_well_index_output(self):
        try:
            from tabulate import tabulate
            print (tabulate([[i, f"{self.name_locus_mapping[i][0]}{self.name_locus_mapping[i][1]+1}"] for i in range(1, len(self.name_locus_mapping)+1)], headers = ["sample", "well"], tablefmt='fancy_grid'))
        except:
            print ("tabulate module not found. please install tabulate module to show the mapping of the sample index to the well index.")
    
    # median value output
    # usage: median_96well(experiment, channel = "Alexa 647-A", interval = "\t")
    # this function will output the median values of the channel in the 96-well plate format.
    # the output will be saved to the working directory with the name "median_values.txt".
    # the interval between the values can be set by the parameter "interval".
    # the median values will also be returned as a 2D list.
    # this function only works for 96-well plate.
    def median_96well(self, experiment, channel = "Alexa 647-A", interval = "\t"):
        median_values_raw = experiment.data.groupby('sample')[channel].median()
        median_values = [[0 for i in range(12)] for j in range(8)]
        for key in median_values_raw.index:
            median_values[ord(self.name_locus_mapping[key][0])-ord('A')][self.name_locus_mapping[key][1]] = median_values_raw[key]
        file = open(self.working_dir + "median_values.txt", "w")
        for i in median_values:
            for j in i:
                file.write(str(j) + interval)
            file.write("\n")
        file.close()
        print (f"median values of {channel} output to {self.working_dir}median_values.txt")
        return median_values

# gate out outliers
def gate_outliers(experiment, channel = "Alexa 647-A", if_plot = False, low = 2000, high = 2147483647):
    gate = cytoflow.RangeOp(name = "de_outliers", channel = channel, low = low, high = high)
    if if_plot:
        cytoflow.HistogramView(channel = channel, scale = "log", huefacet = "well").plot(experiment, num_bins = 200, density = True, title = f"Histogram of {channel}")
    return gate.apply(experiment).subset("de_outliers", True)

# assist auto gate
assistants = []
def assist_gate(experiment, channel = "FSC 488/10-A", if_plot = False, low = 2000, high = 2147483647):
    global assistants
    gate = cytoflow.RangeOp(name = f"assist{len(assistants)+1}", channel = channel, low = low, high = high)
    assistants.append(gate)
    if if_plot:
        gate.default_view(scale = "log").plot(experiment, title = channel)
    return gate.apply(experiment).subset(f"assist{len(assistants)}", True)  # without +1 because already conducted appending

# plot FSC-A ~ SSC-A
def FSC_SSC_ploting(experiment, type = "scatter"):
    if type == "scatter":
        cytoflow.ScatterplotView(xchannel = "FSC 488/10-A", ychannel = "SSC 488/10-A", yscale = "linear", xscale = "linear", huefacet = 'sample').plot(experiment, alpha = 0.005)
    elif type == "density":
        cytoflow.DensityView(xchannel = "FSC 488/10-A", ychannel = "SSC 488/10-A", yscale = "linear", xscale = "linear").plot(experiment)
    return

# gate FSC-A ~ SSC-A
def gate_FSC_SSC(experiment, vertices, if_plot = False):
    gate = cytoflow.PolygonOp(name = "FSCA_SSCA", xchannel = "FSC 488/10-A", ychannel = "SSC 488/10-A", vertices = vertices)
    gatedDATA = gate.apply(experiment)
    if if_plot:
        gate.default_view(xscale = "linear", yscale = "linear", huefacet = 'sample').plot(gatedDATA, s = 10, alpha = 0.1)
    return gatedDATA.subset("FSCA_SSCA", True)

# auto gate FSC-A ~ SSC-A
def auto_gate_FSC_SSC(experiment, if_plot = False, keep = 0.6):
    dens_gate = cytoflow.DensityGateOp(name = "auto_FSCA_SSCA", xchannel = "FSC 488/10-A", ychannel = "SSC 488/10-A", xscale = "linear", yscale = "linear", keep = keep)
    dens_gate.estimate(experiment)
    if if_plot:
        dens_gate.default_view().plot(experiment)
    return dens_gate.apply(experiment).subset("auto_FSCA_SSCA", True)

# plot FSC-A ~ FSC-H
def FSC_A_H_ploting(experiment):
    cytoflow.ScatterplotView(xchannel = "FSC 488/10-A", ychannel = "FSC 488/10-H", yscale = "linear", xscale = "linear", huefacet = 'sample').plot(experiment, alpha = 0.005)
    return

# gate FSC-A ~ FSC-H
def gate_FSC_A_H(experiment, vertices, if_plot = False):
    gate = cytoflow.PolygonOp(name = "FSCA_FSCH", xchannel = "FSC 488/10-A", ychannel = "FSC 488/10-H", vertices = vertices)
    gatedDATA = gate.apply(experiment)
    if if_plot:
        gate.default_view(xscale = "linear", yscale = "linear", huefacet = 'sample').plot(gatedDATA, s = 10, alpha = 0.1)
    return gatedDATA.subset("FSCA_FSCH", True)

# auto gate FSC-A ~ FSC-H
def auto_gate_FSC_A_H(experiment, if_plot = False, keep = 0.9):
    dens_gate = cytoflow.DensityGateOp(name = "auto_FSCA_FSCH", xchannel = "FSC 488/10-A", ychannel = "FSC 488/10-H", xscale = "linear", yscale = "linear", keep = keep)
    dens_gate.estimate(experiment)
    if if_plot:
        dens_gate.default_view().plot(experiment)
    return dens_gate.apply(experiment).subset("auto_FSCA_FSCH", True)

# make scatter plot
def log_density_plot(experiment, xchannel = "FSC 488/10-A", ychannel = "FSC 488/10-H", huefacet = "", **kwargs):
    if huefacet != "":
        cytoflow.DensityView(xchannel = xchannel, ychannel = ychannel, xscale = "log", yscale = "log", huefacet = huefacet).plot(experiment, **kwargs)
    else:
        cytoflow.DensityView(xchannel = xchannel, ychannel = ychannel, xscale = "log", yscale = "log").plot(experiment, **kwargs)
    return

# plot histogram of a channel
def channel_histogram(experiment, channel = "Alexa 647-A", huefacet = "sample", bins = 200, title = "", **kwargs):
    if title == "":
        if huefacet == "conc":
            cytoflow.HistogramView(channel = channel, scale = "log", huefacet = huefacet, huescale = "log").plot(experiment, num_bins = bins, alpha = 0.9, density = True, **kwargs)
        else:
            cytoflow.HistogramView(channel = channel, scale = "log", huefacet = huefacet).plot(experiment, num_bins = bins, density = True, **kwargs)
        return
    else:
        if huefacet == "conc":
            cytoflow.HistogramView(channel = channel, scale = "log", huefacet = huefacet, huescale = "log").plot(experiment, num_bins = bins, alpha = 0.9, density = True, title = title, **kwargs)
        else:
            cytoflow.HistogramView(channel = channel, scale = "log", huefacet = huefacet).plot(experiment, num_bins = bins, density = True, title = title, **kwargs)
        return

# get sub set by char
def subset_by_char(exper, char):
    subset = exper.subset("char", char)
    return subset

# get sub set by num
def subset_by_num(exper, num):
    subset = exper.subset('num', num)
    return subset

# get sub set by well index
def subset_by_well(exper, well):
    if (type(well) == str):
        subset = exper.subset('well', well)
        return subset
    elif (type(well) == list):
        subset = exper.clone(deep = False)
        subset.data = exper.data[exper.data['well'].isin(well)]
        subset.data.reset_index(drop = True, inplace = True)
        return subset
    else:
        raise TypeError("well should be a string or a list of strings, but got " + str(type(well)))

# get count of cells in the experiment
def cell_count(experiment):
    return len(experiment.data)

# 2D Density plot
def density_plot(experiment, xchannel = "FSC 488/10-A", ychannel = "FSC 488/10-H", xscale = "log", yscale = "log", cmap = 'magma_r', **kwargs):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='scatter_density')
    DensityPlot = ax.scatter_density(experiment.data[xchannel], experiment.data[ychannel], cmap = cmap, **kwargs)
    ax.set_box_aspect(1)
    ax.set_xlabel(xchannel)
    ax.set_ylabel(ychannel)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)

    fig.colorbar(DensityPlot, ax=ax, label='Density')
    ax.grid(False)
    return fig, ax

# Draggable Polygon for interactive polygon gate creation
from matplotlib.widgets import Button
class DraggablePolygon:
    def __init__(self, fig, ax):
        self.fig = fig
        self.ax = ax
        plt.subplots_adjust(bottom=0.2)  # Make space for button
        self.points = []
        self.dragging_point = None
        self.epsilon = 10  # pixel distance to select a point
        self.poly = None
        self.scatter = None
        self.interactive_enabled = False

        # Add button for switching interaction
        self.ax_button = plt.axes([0.4, 0.05, 0.2, 0.075])
        self.button = Button(self.ax_button, 'Enable Interaction')
        self.button.on_clicked(self.toggle_interaction)

        self.cid_click = self.fig.canvas.mpl_connect('button_press_event', self.on_click)
        self.cid_release = self.fig.canvas.mpl_connect('button_release_event', self.on_release)
        self.cid_motion = self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.cid_close = self.fig.canvas.mpl_connect('close_event', self.on_close)

        self.ax.set_title(self.get_title())
        plt.show()

    def get_title(self):
        return ("[ON] " if self.interactive_enabled else "[OFF] ") + "Click to add points. Drag to move. Close window to continue and output points."

    def toggle_interaction(self, event):
        self.interactive_enabled = not self.interactive_enabled
        self.button.label.set_text('Disable Interaction' if self.interactive_enabled else 'Enable Interaction')
        self.ax.set_title(self.get_title())
        self.fig.canvas.draw_idle()

    def on_click(self, event):
        if not self.interactive_enabled or event.inaxes != self.ax:
            return
        # Check if click is near an existing point
        for i, (x, y) in enumerate(self.points):
            if (abs(event.x - self.ax.transData.transform((x, y))[0]) < self.epsilon and
                abs(event.y - self.ax.transData.transform((x, y))[1]) < self.epsilon):
                self.dragging_point = i
                return
        # Otherwise, add a new point
        self.points.append([event.xdata, event.ydata])
        self.update_plot()

    def on_release(self, event):
        self.dragging_point = None

    def on_motion(self, event):
        if not self.interactive_enabled or self.dragging_point is None or event.inaxes != self.ax:
            return
        self.points[self.dragging_point] = [event.xdata, event.ydata]
        self.update_plot()

    def update_plot(self):
        # self.ax.clear()
        self.ax.set_title(self.get_title())
        if self.points:
            if len(self.points) >= 1 and self.scatter is None:
                xs, ys = zip(*self.points)
                self.scatter = self.ax.scatter(xs, ys, color='red', zorder=2)
            elif self.scatter is not None:
                xs, ys = zip(*self.points)
                self.scatter.set_offsets(np.c_[xs, ys])
            if len(self.points) > 1 and self.poly is None:
                poly_points = self.points + [self.points[0]]
                px, py = zip(*poly_points)
                self.poly, = self.ax.plot(px, py, marker='o', color='blue', zorder=1)
            elif self.poly is not None:
                poly_points = self.points + [self.points[0]]
                px, py = zip(*poly_points)
                self.poly.set_data(px, py)
        self.fig.canvas.draw()

    def on_close(self, event):
        print("Polygon points:")
        print(self.points)

# polygon gate
polygon_gate_number = 0
def polygon_gate(experiment, xchannel = "FSC 488/10-A", ychannel = "SSC 488/10-A", xscale = "linear", yscale = "linear", cmap='magma_r', vertices = [], if_plot = False):
    # construct a polygon gate if vertices is provided
    global polygon_gate_number
    if vertices != []:
        gate = cytoflow.PolygonOp(name = f"polygon_gate_{polygon_gate_number}", xchannel = xchannel, ychannel = ychannel, vertices = vertices)
        gatedDATA = gate.apply(experiment)
        polygon_gate_number += 1
        if if_plot:
            gate.default_view(xscale = xscale, yscale = yscale).plot(gatedDATA, s = 10, alpha = 0.1)
        return gatedDATA.subset(f"polygon_gate_{polygon_gate_number-1}", True)
    else:
        fig, ax = density_plot(experiment, xchannel = xchannel, ychannel = ychannel, xscale = xscale, yscale = yscale, cmap=cmap)
        draggable_polygon = DraggablePolygon(fig, ax)
        gate = cytoflow.PolygonOp(name = f"polygon_gate_{polygon_gate_number}", xchannel = xchannel, ychannel = ychannel, vertices = draggable_polygon.points)
        gatedDATA = gate.apply(experiment)
        polygon_gate_number += 1
        if if_plot:
            gate.default_view(xscale = xscale, yscale = yscale).plot(gatedDATA, s = 10, alpha = 0.1)
        return gatedDATA.subset(f"polygon_gate_{polygon_gate_number-1}", True)
