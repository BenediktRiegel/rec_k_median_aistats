#!/usr/bin/python3

import matplotlib.pyplot as plt
import numpy as np
import sys
import glob
import pathlib

path = str(pathlib.Path(__file__).parent.absolute()).replace('/..', '')

path_Rec = path + '/EvaluationResults/RecLocalSearch/'
path_kUFL = path + '/EvaluationResults/kUFLLocalSearch/'


def get_k_lam(files):
    ks = set()
    lams = set()
    for file in files:
        file = file.split('_')
        for el in file:
            el = el.split('=')
            if el[0] == 'lam':
                lams.add(float(el[1]))
            elif el[0] == 'k':
                ks.add(int(el[1]))
    return ks, lams


def min_f_str(f):
    f = str(f).split('.')
    all_zero = True
    if len(f) == 1:
        return '.'.join(f)
    for el in f[1]:
        if el != '0':
            all_zero = False
    if all_zero:
        return f[0]
    return '.'.join(f)


def get_content(data_points, path, ks, lams):
    content = {k: {lam: {data_point: [] for data_point in data_points} for lam in lams} for k in ks}
    for k in ks:
        for lam in lams:
            f = open(path + 'k=' + str(k) + '_lam=' + min_f_str(lam) + '.txt')
            for line in f.readlines():
                line = line.split('_')
                for entry in line:
                    entry = entry.split('=')
                    if entry[0] in data_points:
                        if entry[0] == 'serviceCost' or entry[0] == 'otherCost':
                            content[k][lam][entry[0]].append(float(entry[1]))
                        elif entry[0] == 'S':
                            S = entry[1][1:-1].split(',')
                            S = [int(el) for el in S]
                            S.sort()
                            content[k][lam][entry[0]].append(S)
                        elif entry[0] == 'duration':
                            content[k][lam][entry[0]].append(float(entry[1].replace('ms', '')))
                        else:
                            content[k][lam][entry[0]].append(entry[1])
    return content


def read_data(data_points):
    rec_files = [file_path.replace(path_Rec, '').replace('.txt', '') for file_path in glob.glob(path_Rec + '*.txt')]
    (ks, lams) = get_k_lam(rec_files)
    rec_content = get_content(data_points, path_Rec, ks, lams)
    kufl_content = get_content(data_points, path_kUFL, ks, lams)
    return list(ks), list(lams), rec_content, kufl_content


def meandiff_plot():
    data_points = ['serviceCost', 'otherCost']
    (ks, lams, rec_content, kufl_content) = read_data(data_points)
    lams.sort()
    ks.sort()
    for k in ks:
        meandiff = []
        for lam in lams:
            rec_mean = np.add(np.array(rec_content[k][lam]['serviceCost']), np.array(rec_content[k][lam]['otherCost']))
            rec_mean = rec_mean.mean()
            kufl_mean = np.add(np.array(kufl_content[k][lam]['serviceCost']), np.array(kufl_content[k][lam]['otherCost']))
            kufl_mean = kufl_mean.mean()
            meandiff.append(rec_mean - kufl_mean)
        plt.plot(lams, meandiff, label='meandiff', marker='o')
    # plt.title('Mean difference Rec - kUFL')
    plt.xlabel('λ')
    plt.ylabel('Cost(S)')
    plt.legend()
    plt.xticks(lams)


def minmaxdiff_plot():
    data_points = ['serviceCost', 'otherCost']
    (ks, lams, rec_content, kufl_content) = read_data(data_points)
    lams.sort()
    ks.sort()
    for k in ks:
        meandiff = []
        for lam in lams:
            rec_min = np.add(np.array(rec_content[k][lam]['serviceCost']), np.array(rec_content[k][lam]['otherCost']))
            rec_min = rec_min.min()
            kufl_max = np.add(np.array(kufl_content[k][lam]['serviceCost']),
                               np.array(kufl_content[k][lam]['otherCost']))
            kufl_max = kufl_max.max()
            if k == 8 and lam == 1:
                print(np.array(kufl_content[k][lam]['serviceCost']))
                print(np.array(kufl_content[k][lam]['otherCost']))
                print(np.add(np.array(kufl_content[k][lam]['serviceCost']),
                               np.array(kufl_content[k][lam]['otherCost'])))
                print(kufl_max)
            meandiff.append(rec_min - kufl_max)
        plt.plot(lams, meandiff, label='minmax', marker='o')
    # plt.title('difference min Rec - max kUFL')
    plt.xlabel('λ')
    plt.ylabel('Cost(S)')
    plt.legend()
    plt.xticks(lams)


def minmindiff_plot():
    data_points = ['serviceCost', 'otherCost']
    (ks, lams, rec_content, kufl_content) = read_data(data_points)
    lams.sort()
    ks.sort()
    for k in ks:
        meandiff = []
        for lam in lams:
            rec_min = np.add(np.array(rec_content[k][lam]['serviceCost']), np.array(rec_content[k][lam]['otherCost']))
            rec_min = rec_min.min()
            kufl_min = np.add(np.array(kufl_content[k][lam]['serviceCost']),
                               np.array(kufl_content[k][lam]['otherCost']))
            kufl_min = kufl_min.min()
            if k == 8 and lam == 1:
                print(np.array(kufl_content[k][lam]['serviceCost']))
                print(np.array(kufl_content[k][lam]['otherCost']))
                print(np.add(np.array(kufl_content[k][lam]['serviceCost']),
                               np.array(kufl_content[k][lam]['otherCost'])))
                print(kufl_min)
            meandiff.append(rec_min - kufl_min)
        plt.plot(lams, meandiff, label='minmindiff', marker='o')
    # plt.title('difference min Rec - min kUFL')
    plt.xlabel('λ')
    plt.ylabel('Cost(S)')
    plt.legend()
    plt.xticks(lams)


def mean_plot():
    data_points = ['serviceCost', 'otherCost']
    (ks, lams, rec_content, kufl_content) = read_data(data_points)
    lams.sort()
    ks.sort()
    for k in ks:
        rec_means = []
        rec_stds = []
        kufl_means = []
        kufl_stds = []
        for lam in lams:
            rec_mean = np.add(np.array(rec_content[k][lam]['serviceCost']), np.array(rec_content[k][lam]['otherCost']))
            rec_means.append(rec_mean.mean())
            rec_stds.append(rec_mean.std())
            kufl_mean = np.add(np.array(kufl_content[k][lam]['serviceCost']),
                              np.array(kufl_content[k][lam]['otherCost']))
            kufl_means.append(kufl_mean.mean())
            kufl_stds.append(kufl_mean.std())
        plt.errorbar(lams, rec_means, yerr=rec_stds, capsize=5.0, capthick=1.0, label='LocalSearch', color='green', marker='o')
        plt.errorbar(lams, kufl_means, yerr=kufl_stds, capsize=5.0, capthick=1.0, label='kUFL+LS', color='blue', marker='o')
    # plt.title('mean Rec and mean kUFL')
    plt.xlabel('λ')
    plt.ylabel('Cost(S)')
    plt.legend()
    plt.xticks(lams)
    

def min_plot():
    data_points = ['serviceCost', 'otherCost']
    (ks, lams, rec_content, kufl_content) = read_data(data_points)
    lams.sort()
    ks.sort()
    for k in ks:
        rec_mins = []
        kufl_mins = []
        for lam in lams:
            rec_min = np.add(np.array(rec_content[k][lam]['serviceCost']), np.array(rec_content[k][lam]['otherCost']))
            rec_min = rec_min.min()
            rec_mins.append(rec_min)
            kufl_min = np.add(np.array(kufl_content[k][lam]['serviceCost']),
                              np.array(kufl_content[k][lam]['otherCost']))
            kufl_min = kufl_min.min()
            kufl_mins.append(kufl_min)
        plt.plot(lams, rec_mins, label='LocalSearch', color='green', linestyle='--', marker='o')
        plt.plot(lams, kufl_mins, label='kUFL+LS', color='blue', linestyle='--', marker='o')
    # plt.title('min Rec and min kUFL')
    plt.xlabel('λ')
    plt.ylabel('Cost(S)')
    plt.legend()
    plt.xticks(lams)


def meanratio_plot():
    data_points = ['serviceCost', 'otherCost']
    (ks, lams, rec_content, kufl_content) = read_data(data_points)
    lams.sort()
    ks.sort()
    for k in ks:
        meanratio = []
        for lam in lams:
            rec_mean = np.add(np.array(rec_content[k][lam]['serviceCost']), np.array(rec_content[k][lam]['otherCost']))
            rec_mean = rec_mean.mean()
            kufl_mean = np.add(np.array(kufl_content[k][lam]['serviceCost']),
                               np.array(kufl_content[k][lam]['otherCost']))
            kufl_mean = kufl_mean.mean()
            meanratio.append(rec_mean / kufl_mean)
        plt.plot(lams, meanratio, label='meanratio', marker='o')
    # plt.title('Mean ratio Rec / kUFL')
    plt.xlabel('λ')
    plt.ylabel('Cost Ratio')
    plt.legend()
    plt.xticks(lams)


def separatemeanratio_plot():
    data_points = ['serviceCost', 'otherCost']
    (ks, lams, rec_content, kufl_content) = read_data(data_points)
    lams.sort()
    ks.sort()
    for k in ks:
        meanotherratio = []
        meanserviceratio = []
        for lam in lams:
            other_rec_mean = np.array(rec_content[k][lam]['otherCost'])
            other_rec_mean = other_rec_mean.mean()
            service_kufl_mean = np.array(kufl_content[k][lam]['otherCost'])
            service_kufl_mean = service_kufl_mean.mean()
            meanotherratio.append(other_rec_mean / service_kufl_mean)

            service_rec_mean = np.array(rec_content[k][lam]['serviceCost'])
            service_rec_mean = service_rec_mean.mean()
            service_kufl_mean = np.array(kufl_content[k][lam]['serviceCost'])
            service_kufl_mean = service_kufl_mean.mean()
            meanserviceratio.append(service_rec_mean / service_kufl_mean)
        plt.plot(lams, meanotherratio, label='Rec. Cost', marker='o')
        plt.plot(lams, meanserviceratio, label='Service Cost', marker='o')
        plt.plot(lams, [1 for el in meanotherratio], label='1')
    # plt.title('Separate Mean ratio Rec / kUFL')
    plt.xlabel('λ')
    plt.ylabel('Cost Ratio')
    plt.legend()
    plt.xticks(lams)
    # ticks = [0]
    # for el in lams:
    #     if el >= 1:
    #         ticks.append(el)
    # plt.xticks(ticks)


def count_good_s(A, S):
    sum = 0
    for s in S:
        if np.array_equal(s, A):
            sum += 1
    return sum


def probability_plot():
    data_points = ['S']
    (ks, lams, rec_content, kufl_content) = read_data(data_points)
    lams.sort()
    ks.sort()
    for k in ks:
        rec_probs = []
        kufl_probs = []
        A = np.array(list(range(k)))  # A is the perfect set
        for lam in lams:
            rec_s = np.array(rec_content[k][lam]['S'])
            kufl_s = np.array(kufl_content[k][lam]['S'])
            rec_probs.append(count_good_s(A, rec_s) / float(len(rec_s)))
            kufl_probs.append(count_good_s(A, kufl_s) / float(len(kufl_s)))
        plt.plot(lams, rec_probs, label='LocalSearch', marker='o')
        plt.plot(lams, kufl_probs, label='kUFL+LS', marker='o')
    # plt.title('The probability of ending with a good set')
    plt.xlabel('Sample Amount')
    plt.ylabel('Estimated Pr$\\left(S=O_{RC}\\right)$')
    plt.legend()
    plt.xticks(lams)


def mean_duration_plot(factor=1.):
    data_points = ['duration']
    (ks, lams, rec_content, kufl_content) = read_data(data_points)
    lams.sort()
    ks.sort()
    for k in ks:
        rec_durations = []
        rec_stds = []
        kufl_durations = []
        kufl_stds = []
        for lam in lams:
            rec_duration = np.array(rec_content[k][lam]['duration'])
            kufl_duration = np.array(kufl_content[k][lam]['duration'])
            rec_durations.append(rec_duration.mean() / factor)
            kufl_durations.append(kufl_duration.mean() / factor)
            rec_stds.append(rec_duration.std() / factor)
            kufl_stds.append(kufl_duration.std() / factor)
        plt.errorbar(lams, rec_durations, yerr=rec_stds, capsize=5.0, capthick=1.0, label='LocalSearch', marker='o')
        plt.errorbar(lams, kufl_durations, yerr=kufl_stds, capsize=5.0, capthick=1.0, label='kUFL+LS', marker='o')
    # plt.title('Mean Runtime')
    plt.xlabel('λ')
    if factor == 1000.:
        plt.ylabel('Runtime [s]')
    else:
        plt.ylabel('Runtime [ms]')
    plt.legend()
    plt.xticks(lams)


def mean_duration_diff_plot(factor=1.):
    data_points = ['duration']
    (ks, lams, rec_content, kufl_content) = read_data(data_points)
    lams.sort()
    ks.sort()
    for k in ks:
        duration_diffs = []
        for lam in lams:
            rec_duration = np.array(rec_content[k][lam]['duration'])
            kufl_duration = np.array(kufl_content[k][lam]['duration'])
            duration_diffs.append((kufl_duration.mean() - rec_duration.mean()) / factor)
        plt.plot(lams, duration_diffs, label='kufl_duration - rec_duration=' + str(k), marker='o')
    # plt.title('Difference of both Mean Runtimes')
    plt.xlabel('λ')
    if factor == 1000.:
        plt.ylabel('Runtime [s]')
    else:
        plt.ylabel('Runtime [ms]')
    plt.legend()
    plt.xticks(lams)


def set_fontsize(text=10, title=10, labels=10, xtick=10, ytick=10, legend=10):
    plt.rc('font', size=text)  # controls default text size
    plt.rc('axes', titlesize=title)  # fontsize of the title
    plt.rc('axes', labelsize=labels)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=xtick)  # fontsize of the x tick labels
    plt.rc('ytick', labelsize=ytick)  # fontsize of the y tick labels
    plt.rc('legend', fontsize=legend)  # fontsize of the legend


def plot(ptypes):
    one_plot = False
    set_fontsize(text=22, title=0, labels=22, xtick=17, ytick=17, legend=17)

    for ptype in ptypes:
        if ptype == 'oneplot':
            one_plot = True

    for ptype in ptypes:
        if ptype == 'meandiff':
            meandiff_plot()
        elif ptype == 'minmaxdiff':
            minmaxdiff_plot()
        elif ptype == 'minmindiff':
            minmindiff_plot()
        elif ptype == 'min':
            min_plot()
        elif ptype == 'mean':
            mean_plot()
        elif ptype == 'meanratio':
            meanratio_plot()
        elif ptype == 'separatemeanratio':
            separatemeanratio_plot()
        elif ptype == 'probability':
            probability_plot()
        elif ptype == 'meanduration':
            mean_duration_plot()
        elif ptype == 'meanduration-s':
            mean_duration_plot(factor=1000.)
        elif ptype == 'meandurationdiff':
            mean_duration_diff_plot()
        elif ptype == 'meandurationdiff-s':
            mean_duration_diff_plot(factor=1000.)
        if not one_plot:
            # plt.xlabel('Number of Clients')
            # plt.xlabel('Sample Size')
            ticks = [0, 16, 32, 64, 128]
            plt.xticks(ticks)
            plt.xlabel('λ')
            plt.tight_layout()
            plt.savefig(ptype + '.pdf')
            plt.clf()
            plt.cla()
    if one_plot:
        # plt.title('oneplot')
        plt.savefig('_'.join(ptypes) + '.pdf')


def ready_argv(argv):
    global path_Rec, path_kUFL
    argv = sys.argv[1:]
    if 'random' in argv:
        path_Rec = path + '/EvaluationResults/RecLocalSearchRandomSampling/'
        path_kUFL = path + '/EvaluationResults/kUFLLocalSearchRandomSampling/'
        argv.remove('random')
    return argv


if __name__ == '__main__':
    plot(ready_argv(sys.argv))

