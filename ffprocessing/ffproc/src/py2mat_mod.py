"""
Convenient method to convert Python objects into MATLAB format
This is mainly because MATLAB plotting tools and processing is much
easier and BETTER! :D ... we will see
"""

import scipy.io as sio

def py2mat(pyobject, name='pyobj', filename='out'):
    """
    input:
    pyobject: the Python object that you want to convert to MATLAB
    name: given name (this will be used in MATLAB)
    filename: the name of the file.mat
    :param pyobject:
    :param name:
    :param filename:
    :return:
    """

    sio.savemat(filename, {name: pyobject})
    print "\n\n===================="
    print "converted to %s (MATLAB) and saved in %s.mat" \
                    %(name, filename)
    print "\nFrom matlab consule type:"
    print "load %s.mat" %(filename)
    print "%s" %(name)
    print "===================="
