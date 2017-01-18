#!/usr/bin/env python
# python 2.7


def save_pydata(data, outfile):
    # save py data in pickle format
    # USAGE: save_pydata(python_object, "my_file.pickle") 
    import pickle
    with open(outfile, 'wb') as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)
        print 'Object saved to file:\n', outfile


def load_pydata(infile):
    # open py pickle data
    import pickle
    with open(infile, 'rb') as f:
        # The protocol version used is detected automatically, so we do not
        # have to specify it.
        data = pickle.load(f)
    return data

def my_debugger(vars):
    # starts interactive Python terminal at location in script
    # call with python_functions.my_debugger(globals().copy()) anywhere in your script
    import readline # optional, will allow Up/Down/History in the console
    import code
    # vars = globals().copy() # in python "global" variables are actually module-level
    vars.update(locals())
    shell = code.InteractiveConsole(vars)
    shell.interact()
