#------------------------------------------------------------------------------
# Initial stuff; use functions for 
#------------------------------------------------------------------------------
import numpy as np
import os
import IPython as ipy
# __all__ = [
#     'dict', # dictionary with dot-access for attributes; simple
#     'notebook', # special IPython notebook stuff
#     'flatten', 'unflatten', 'permute', 'unpermute',
#     'match', 'arange', 'range',
#     'slope', 'rolling', 'runslope', 'runmean',
#     ]

#------------------------------------------------------------------------------
# Jupyter Notebook setup
#------------------------------------------------------------------------------
import IPython as ipy
def notebook(stem=None):
    """
    Setup up magic commands and other notebook properties.
    Input...
        stem: the directory stem from HOME directory.
    """
    # Directory
    if stem is None:
        os.chdir(os.environ['HOME'] + '/Desktop')
    else:
        os.chdir(os.environ['HOME'] + '/' + stem.lstrip('/'))

    # CSS modifications; use here, or place in ~/.jupyter/custom/custom.css
    # import IPython.core.display as display # use these to edit CSS
    # ipy.core.display.display(ipy.core.display.HTML("<style>.container { width:100% !important; }</style>"))
    # display.display(display.HTML("<style>"
    #     + "#notebook { padding-top:0px !important; }" # hashtags refer to ids; dots refer to classes
    #     + ".container { width:100% !important; } "
    #     + ".end_space { min-height:0px !important; }"
    #     + "</style>"))
    # # icolor, ncolor = "#2b2b25"
    # # icolor, ncolor = "#f2f2f2", "#f2f2f2"
    # icolor, ncolor = "#f6f6f6", "#f6f6f6"
    # display.display(display.HTML("<style>"
    #     + ".edit_mode .cell.selected .CodeMirror-focused.cm-fat-cursor { background-color: " + ncolor + " !important; } " # normal mode
    #     + ".edit_mode .cell.selected .CodeMirror-focused:not(.cm-fat-cursor) { background-color: " + icolor + " !important; }" # edit mode
    #         # the above matches color of ITerm2 settings
    #     + "</style>")) # <style> denotes CSS code block

    # Magic commands
    running = ipy.get_ipython()
    # ipy.core.display.display(ipy.core.display.HTML(html_format))
    running.magic('autosave 120') # autosave every 120 seconds
    running.magic('load_ext autoreload')
    running.magic('autoreload 2')
    running.magic('config InlineBackend.figure_format=\'svg\'')
        # retina probably more space efficient
    # running.magic('config InlineBackend.figure_format=\'svg\'')
    # running.magic('config InlineBackend.print_figure_kwargs=dict(bbox_inches=\'tight\', pad_inches=0.05)')
    running.magic('config InlineBackend.print_figure_kwargs=dict(bbox_inches=None)') #bbox_inches=\'tight\', pad_inches=0.1)')
    running.magic('matplotlib inline')

#------------------------------------------------------------------------------
# General definitions
#------------------------------------------------------------------------------
class Dict(dict):
    """
    Dot notation access to dictionary attributes; EXTREMELY useful, and no
    slower than normal dictionaries.
    """
    # Very simple
    # __getattr__ = dict.get
    # __setattr__ = dict.__setitem__
    # __delattr__ = dict.__delitem__
    # More complex implimentation
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for arg in args:
            if isinstance(arg, dict):
                for k, v in arg.items():
                    self[k] = v
        if kwargs: # if non-empty
            for k, v in kwargs.items():
                self[k] = v

    def __getattr__(self, attr):
        return self.get(attr)

    def __setattr__(self, key, value):
        self.__setitem__(key, value)

    def __setitem__(self, key, value):
        super().__setitem__(key, value)
        self.__dict__.update({key: value})

    def __delattr__(self, item):
        self.__delitem__(item)

    def __delitem__(self, key):
        super().__delitem__(key)
        del self.__dict__[key]

def year(dt):
    """
    Gets year from numpy datetime object.
    """
    return dt.astype('datetime64[Y]').astype(np.int32)+1970 # the astype(int) is actually super fast (ns)
        # and the above in general is faster than list comprehension with container of datetime objects
        # UNIX time starts at 1970

def month(dt):
    """
    Gets month from numpy datetime object.
    """
    return dt.astype('datetime64[M]').astype(np.int32)%12 + 1
        # will convert datetime64 units from [ns] (default) to months, then spit out
        # numerical months relative to year

def arange(min_, *args):
    """
    Duplicate behavior of np.arange, except with inclusive endpoints; dtype is
    controlled very carefully, so should be 'most precise' among min/max/step args.
    If any of min/max/step are floating point, use np.nextafter to increment as small as possible 
    for the np.arange endpoint.
    Input...
        stop
        start, stop, [step]
        ...just like np.arange
    Output...
        the array sequence
    """
    # Optional arguments just like np.arange
    if len(args)==0:
        max_ = min_
        min_ = 0 # this re-assignes the NAME "min_" to 0
        step = 1
    elif len(args)==1:
        max_ = args[0]
        step = 1
    elif len(args)==2:
        max_ = args[0]
        step = args[1]
    else:
        raise ValueError('Function takes from one to three arguments.')
    # All input is integer? Get new "max"
    if min_//1==min_ and max_//1==max_ and step//1==step:
        min_, max_, step = np.int64(min_), np.int64(max_), np.int64(step)
        max_ += 1
    # Input is float or mixed; cast all to float64, then get new "max"
    else:
        min_, max_, step = np.float64(min_), np.float64(max_), np.float64(step)
        max_ += step/2
        # max_ = np.nextafter(max_, np.finfo(np.dtype(np.float64)).max)
            # gives the next FLOATING POINT, in direction of the second argument
            # ...forget this; round-off errors from continually adding step to min mess this up
    return np.arange(min_, max_, step)

def match(v1, v2):
    """
    Match two 1D vectors; will return starting/ending indices, and 
    numpy array of their overlap.
    """
    v1, v2 = np.array(v1), np.array(v2)
    if not np.all(v1==np.sort(v1)) or not np.all(v2==np.sort(v2)):
        raise ValueError('Vectors must be sorted.')
    # Get common minima/maxima
    min12, max12 = max(v1.min(), v2.min()), min(v1.max(), v2.max())
    try:
        min1f, min2f = np.where(v1==min12)[0][0], np.where(v2==min12)[0][0]
        max1f, max2f = np.where(v1==max12)[0][0], np.where(v2==max12)[0][0]
    except IndexError:
        raise ValueError('Vectors do not have matching maxima/minima.')
    slice1, slice2 = slice(min1f, max1f+1), slice(min2f, max2f+1)
    if v1[slice1].size != v2[slice2].size:
        raise ValueError('Vectors are not identical between matching minima/maxima.')
    elif not (v1[slice1]==v2[slice2]).all():
        raise ValueError('Vectors are not identical between matching minima/maxima.')
    return slice1, slice2, v1[slice1]

