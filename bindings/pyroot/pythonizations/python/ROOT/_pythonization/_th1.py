# Author: Enric Tejedor CERN  02/2019

################################################################################
# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

from . import pythonization


def __init__(self, *binnings, counts=None):
    if len(binnings) == 0:
        raise TypeError("Histogram must not have zero dimensions")
    if any(not isinstance(binning, Binning) for binning in binnings):
        raise TypeError("non-keyword arguments to Histogram must all be binnings")
    self.binnings = binnings

    if counts is None:
        counts = numpy.zeros([len(binning) for binning in self.binnings], int)
    if isinstance(counts, numbers.Integral):
        counts = numpy.full([len(binning) for binning in self.binnings], counts)
    if not isinstance(counts, numpy.ndarray):
        counts = numpy.array(counts)
    self.counts = counts

    assert len(binnings) == len(counts.shape)
    assert all(len(binning) == length for binning, length in zip(binnings, counts.shape))

def __repr__(self):
    return "Histogram({0})".format(", ".join(repr(binning) for binning in self.binnings))

def __getitem__(self, slices):
    if not isinstance(slices, tuple):
        slices = (slices,)
    if len(slices) > len(self.binnings):
        raise IndexError("too many slices for Histogram of dimension {0}".format(len(self.binnings)))
    slices = slices + (slice(None),) * (len(self.binnings) - len(slices))

    binnings = []
    counts = self.counts
    for item, binning in zip(slices, self.binnings):
        if not isinstance(item, slice):
            raise TypeError("only slices are allowed")
        start, stop, step = item.start, item.stop, item.step

        # the current axis depends on the dimensionality of the output
        axis = len(binnings)

        # slice the binning, creating under/overflows if necessary
        if start is not None or stop is not None:
            binning, counts = binning.sliced(start, stop, axis, counts, wantflows=getattr(step, "wantflows", True))

        # apply the project/rebin function
        if step is not None:
            if not callable(step):
                raise TypeError("when slicing a Histogram, the slice's third argument must be callable")
            original_shape = counts.shape

            # signature: (binning for axis=axis, axis number, whole counts array) → new binning (or None), new counts
            binning, counts = step(binning, axis, counts)

            # ensure that the (possibly user-supplied) function is sane
            if binning is None:
                assert original_shape[:axis] + original_shape[axis + 1:] == counts.shape
            else:
                assert original_shape[:axis] == counts.shape[:axis]
                assert original_shape[axis + 1:] == counts.shape[axis + 1:]
                assert len(binning) == counts.shape[axis]

        # also determines whether the axis number will be increased in the next round
        if binning is not None:
            binnings.append(binning)

    return Histogram(*binnings, counts=counts)

def __eq__(self, other):
    return isinstance(other, Histogram) and self.binnings == other.binnings and numpy.array_equal(self.counts, other.counts)

def __ne__(self, other):
    return not self.__eq__(other)

def loc(x):
    """When used in the start or stop of a Histogram's slice, x is taken to be the position in data coordinates."""
    return lambda binning, isleft: binning.index(x, isleft=isleft, clip=True)

def project(binning, axis, counts):
    """When used in the step of a Histogram's slice, project sums over and eliminates what remains of the axis after slicing."""
    return None, numpy.add.reduce(counts, axis=axis)

project.wantflows = False

def rebin(factor):
    """When used in the step of a Histogram's slice, rebin(n) combines bins, scaling their widths by a factor of n. If the number of bins is not divisible by n, the remainder is added to the overflow bin."""
    def impl(binning, axis, counts):
        if isinstance(binning, Regular):
            indexes = (numpy.arange(0, binning.num, factor),)

            num, remainder = divmod(binning.num, factor)
            high, hasover = binning.high, binning.hasover

            if binning.hasunder:
                indexes[0][:] += 1
                indexes = ([0],) + indexes

            if remainder == 0:
                if binning.hasover:
                    indexes = indexes + ([binning.num + int(binning.hasunder)],)
            else:
                high = binning.left(indexes[-1][-1])
                hasover = True

            binning = Regular(num, binning.low, high, hasunder=binning.hasunder, hasover=hasover)
            counts = numpy.add.reduceat(counts, numpy.concatenate(indexes), axis=axis)
            return binning, counts

        else:
            raise NotImplementedError(type(binning))

    return impl

rebin.wantflows = True

# Multiplication by constant

def _imul(self, c):
    # Parameters:
    # - self: histogram
    # - c: constant by which to multiply the histogram
    # Returns:
    # - A multiplied histogram (in place)
    self.Scale(c)
    return self


def FromNumpy(x, w, name, title, nbins, xmin, xmax, size):
    import ROOT
    obj = ROOT.TH1D(name, title, nbins, xmin, xmax)
    obj.FillN(size, x, w)
    return obj


@pythonization('TH1')
def pythonize_th1(klass):
    # Parameters:
    # klass: class to be pythonized
    # Support hist *= scalar
    klass.__imul__ = _imul
    klass.FromNumpy = FromNumpy
