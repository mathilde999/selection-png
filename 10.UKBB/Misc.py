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

def systemcode2output(command):
    """
    If the command are gicen it will crate a fake output file which then later can be easily accessed by python. in case
    of error it will print the error and exit
    :param command: The fucking syustem code. for example ms code etc
    :return: Will return a fake input file (generally a memory location). Which then can be used as iterator of line
    read. for example:
    for line in systemcode2output(command):
        do something
    Good as We dont have to save it in the hd
    """
    import subprocess
    fid = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, universal_newlines=True)
    return fast_reading_by_line_iterator(fid.stdout)

def fast_reading_by_line_iterator(fileObj):
    """
    Fastest way I have found till now to read files line by line
    f = gzip.open(filename,'rb')
    for line in fast_reading_by_line(f):
        DO SOME THING
    f.close()

    :param fileObj: created file object by open and other methods
    :return: will return line by line
    """
    while True:
        data = fileObj.readline().rstrip()
        if not data:
            break
        yield data

def filenamewithoutextension(filepath):
    """
    As the name sugges it will remove the file extension from the full filename. Addtionally if its in a path it will
    remove the path as well. IF it has multiple dot will only remove the first one (i.e. bla.vcf.gz will give bla.vcf)
    :param filepath: The file path
    :return: will give only file name without extension.
    """
    filename = gettingfilename(filepath)
    extension = gettingextension(filepath)
    if len(extension) > 0:
        return filename[:-(len(extension) + 1)]
    else:
        return filename

def gettingfilename(filepath):
    """
        It wil give back the name of the file from the filepath
    :return:    Name of the file
    """
    if "/" in filepath:
        words = filepath.split("/")
        filename = words[len(words) - 1]
        return filename
    else:
        return filepath


def gettingextension(filepath):
    """
    As the name suggest from the filepath it will give the file extension. Remember if the filename do not have dot (.)
    it will give back the whole filename
    :param filepath: The whole file path of the file whose extension we wanted to get
    :return:the extension itself
    """
    import re
    filename = gettingfilename(filepath)
    if re.search(".", filename):
        splitfilename = filename.split(".")
        return splitfilename[len(splitfilename) - 1]
    else:
        return ""