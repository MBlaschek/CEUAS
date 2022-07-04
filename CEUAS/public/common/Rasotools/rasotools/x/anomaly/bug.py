import numpy


# pythran export bug(float[])
def bug(x):

	ints=numpy.zeros(5,numpy.int32)
        ints[0]=numpy.int32(numpy.floor(x))

	return ints

print bug(2.1)

