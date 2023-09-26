import os


def args_valid_file(parser, arg):
    """
    This will check with in the argparse and report if the file exist or not. should be used in the argparse command
    parser.add_argument('file',type=lambda x: Misc.args_valid_file(parser, x))
    Args:
        parser: argparse.ArgumentParser which is the main for argparse command
        arg: the path of the command itself

    Returns: if the file exist it will return the path if not it will print an error message

    """
    if arg:
        if not os.path.exists(arg):
            parser.error("The file %s does not exist!" % arg)
    return arg


def adding_pop(alldata, popfilepath=None, popdf=None):
    """
    This will return a column with pop information on a dataframe
    :param alldata: the data, where you need at least one column with inds information and it should be the first column
    :param popfilepath: the path of pop infor. The file should have first column with inds and the second is for pop.
    No header
    :param popdf: in case already have dataframe loaded for populations. it will ignore the popfile path
    :return: will return a column (Series) of exact size of input (alldata) so that it can be concatenate with the
    data itself
    """
    import numpy
    import pandas
    if popdf is None:
        popfile = pandas.read_csv(popfilepath, usecols=[0, 1], header=None, names=['inds', 'pop'],
                                  delim_whitespace=True, dtype=object)
    else:
        popfile = popdf.iloc[:, :2]
        popfile.columns = ['inds', 'pop']
    alldata = alldata[:]
    alldata = alldata.assign(pop=numpy.nan)
    for index, row in alldata.iterrows():
        try:
            alldata.loc[index, "pop"] = list(popfile[popfile["inds"] == row.iloc[0]]["pop"].values)[0]
        except IndexError:
            pass
    alldata['pop'] = alldata['pop'].replace('nan', numpy.nan)
    return alldata['pop']


def file_existence_checker(names, files):
    """
    As most of the clasess needs different type of check of existence of the file. It will check if the file exist and
    then return appropiate error message to print.
    One way to use the output
    tobeprint, absentfiles = Misc.file_existence_checker(names, files)
    if len(absentfiles) > 0:
        print (tobeprint)
        sys.exit(1)
    :param names: the name or the purpose of the file
    :param files: The file which we are checking
    :return: If exist the files it will return nothing. If any file does not exist it will return the indexes which file
     do not exist with a printable command saying which file (with the name or purpose) does not exist
    """
    import os
    printer = []
    absence = []
    for index, file in enumerate(files):
        if not os.path.isfile(file):
            printerlinetemp = ["The", names[index], "file could not be found:", file]
            printer.append(joinginglistbyspecificstring(printerlinetemp))
            absence.append(file)
    return joinginglistbyspecificstring(printer, "\n"), absence


def joinginglistbyspecificstring(listinput, string=" "):
    """
    It join a list with the given string
    :param listinput: the input list
    :param string: the string with which it has to join. by default is space
    :return: the joined string
    """
    listinput = [x for x in listinput if x is not None]  ##to remove none
    listinput = list(map(str, listinput))
    return string.join(listinput)

