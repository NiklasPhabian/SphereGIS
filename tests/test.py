import sphereGIS 
import numpy
import pandas


class MainTest(unittest.TestCase):
    
    def test1(self):
        df = pandas.read_csv('data/brazil.csv', header=None)
        x = numpy.array(df[0])
        y = numpy.array(df[1])
        z = numpy.array(df[2])

        out = numpy.full(x.shape, [-1], dtype=numpy.int32)


        print(sphereGIS.xyz2convex(x,y,z))
