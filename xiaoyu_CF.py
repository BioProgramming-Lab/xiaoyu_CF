import cytoflow
import glob
import os.path
import numpy as np

import warnings
warnings.filterwarnings("ignore")

class xiaoyu_Expr:

    def __init__(self, working_dir : str):
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

        # classify samples with name (A1, A2, A3, ..., B1, B2, B3, ...)
        file_names = []
        for i in filelist:
            file_names.append(os.path.splitext(os.path.split(i)[-1])[0])
            file_names[-1] = file_names[-1].split(" ")[0]   # filter out the well name. e.g. "A1 EGFR.fcs" -> "A1"
        for i in range(1, len(file_names)+1):
            self.name_locus_mapping[i] = [file_names[i-1][0], int(file_names[i-1][1:])-1]
        #print (f"mapping of sample name to char and num: {name_locus_mapping}")
        # generate experiment instance
        tubes = []
        count = 1
        for i in filelist:
            tube = cytoflow.Tube(file = i, conditions = {"sample" : count, "char" : self.name_locus_mapping[count][0], "num" : str(self.name_locus_mapping[count][1]+1), "well" : self.name_locus_mapping[count][0]+str(self.name_locus_mapping[count][1]+1)})
            tubes.append(tube)
            count += 1
        experiment = cytoflow.ImportOp(conditions = {"sample" : "int", "char" : "category", "num" : "category", "well" : "str"}, tubes = tubes).apply()
        print (f"properties of the experiment: {experiment.data.columns}")
        self.expr = experiment

    # add concetration condition to the experiment
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
    def add_condition_with_sample_array(self, condition_name, dtype, conditions):
        condition_list = []
        for i in self.expr.data["sample"]:
            condition_list.append( conditions[i-1] )
        condition_list = np.array(condition_list)
        self.expr.add_condition(name = condition_name, dtype = dtype, data = condition_list)

    # median value output
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
        return

# gate out outliers
def gate_outliers(experiment, channel = "Alexa 647-A", if_plot = False, low = 2000, high = 2147483647):
    gate = cytoflow.RangeOp(name = "de_outliers", channel = channel, low = low, high = high)
    if if_plot:
        gate.default_view().plot(experiment)
    return gate.apply(experiment).subset("de_outliers", True)

# assist auto gate
def assist_gate(experiment, channel = "FSC 488/10-A", if_plot = False, low = 2000, high = 2147483647):
    gate = cytoflow.RangeOp(name = "assist", channel = channel, low = low, high = high)
    if if_plot:
        gate.default_view().plot(experiment)
    return gate.apply(experiment).subset("assist", True)

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

# plot histogram of a channel
def channel_histogram(experiment, channel = "Alexa 647-A", huefacet = "sample", bins = 200, title = ""):
    if title == "":
        if huefacet == "conc":
            cytoflow.HistogramView(channel = channel, scale = "log", huefacet = huefacet, huescale = "log").plot(experiment, num_bins = bins, alpha = 0.9, density = True)
        else:
            cytoflow.HistogramView(channel = channel, scale = "log", huefacet = huefacet).plot(experiment, num_bins = bins, density = True)
        return
    else:
        if huefacet == "conc":
            cytoflow.HistogramView(channel = channel, scale = "log", huefacet = huefacet, huescale = "log").plot(experiment, num_bins = bins, alpha = 0.9, density = True, title = title)
        else:
            cytoflow.HistogramView(channel = channel, scale = "log", huefacet = huefacet).plot(experiment, num_bins = bins, density = True, title = title)
        return

# get sub set by char
def subset_by_char(exper, char):
    subset = exper.subset("char", char)
    return subset

# get sub set by num
def subset_by_num(exper, num):
    subset = exper.subset('num', num)
    return subset