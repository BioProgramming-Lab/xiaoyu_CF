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