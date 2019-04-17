class NoParameterError(Exception):

    def __init__(self, parameter, prompt):
        self.par = parameter
        self.pro = prompt

    def __str__(self):
        return repr("Value for parameter '%s' (%s) is required!!!" % (self.par, self.pro))


class WrongValueError(Exception):

    def __init__(self, parameter, value, prompt):
        self.par = parameter
        self.val = value
        self.pro = prompt

    def __str__(self):
        return repr("Value %s for parameter '%s' is wrong!!! %s" % (self.val, self.par, self.pro))


class OtherError(Exception):

    def __init__(self, prompt):
        self.pro = prompt

    def __str__(self):
        return repr(self.pro)


class DirectoryError(Exception):

    def __init__(self, parameter):
        self.par = parameter

    def __str__(self):
        return 'After value of parameter %s there should be a directory to folder which contains it.' % self.par


class WrongSubscripts(Exception):

    def __init__(self, oldi, oldr, oldc, i, r, c):
        self.pro = 'Wrong subscripts! Previous %d:<%d, %d>, next %d:<%d, %d>' % (oldi, oldr, oldc, i, r, c)

    def __str__(self):
        return repr(self.pro)


class DiagnoseOverwriting(Exception):

    def __init__(self, dataset, pid):
        self.pro = 'Diagnosis for patient %s from %s has already been written!' % (pid, dataset)

    def __str__(self):
        return repr(self.pro)


class FilesError(Exception):

    def __init__(self, type, i, c):
        if type == 'id':
            self.pro = "Different ID of patient in the %d row for chr %d" % (i, c)
        elif type == 'number':
            self.pro = "Different number of patients: chr = %d, number of patients = %d" % (c, i)
        else:
            self.pro = "There is no such type (%s) of FilesError!" % type

    def __str__(self):
        return repr(self.pro)


class KeyOverwriting(Exception):

    def __init__(self, key):
        self.pro = 'Error: key %s has already been in the dictionary!' % key

    def __str__(self):
        return repr(self.pro)


class PlinkWrongValue(Exception):

    def __init__(self, snp, ch, pat, a, aa):
        self.pro = 'Error: SNP number %s from chr %s for patient %d has unpredicted values - %s and %s!' % \
                   (snp, ch, pat, a, aa)

    def __str__(self):
        return repr(self.pro)


class NoFileError(Exception):

    def __init__(self, file):
        self.pro = "There is no %s file in input directory. It is required for analysis!" % file

    def __str__(self):
        return repr(self.pro)


class SNPReferenceError(Exception):

    def __init__(self, ch, position, ref1, ref2):

        self.pro = 'Different reference value for SNP on position %d from chromosome %d. Set 1: %s, set 2: %s\n' % \
                   (position, ch, ref1, ref2)

    def __str__(self):
        return repr(self.pro)


class FileOverwriteError(Exception):

    def __init__(self, filename):

        self.pro = 'File %s filename already exist!' % filename

    def __str__(self):
        return repr(self.pro)


class NoSNPFound(Exception):

    def __init__(self, ch, perc):
        self.pro = 'Warning: no SNPs were chosen from chromosome %d using Boruta with class perc = %d' % (ch, perc)

    def __str__(self):
        return repr(self.pro)


class WrongParameterName(Exception):

    def __init__(self, parameter):
        self.pro = "There is no command parameter named '%s'!" % parameter

    def __str__(self):
        return repr(self.pro)
