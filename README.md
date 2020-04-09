# SphereGIS


# The goal of SphereGIS:

## Given are:

* A very large set of points p
* A spherical polygon P defined by its nodes N and its edges E

## Goal:
* Find is the subset of p that are within P.
* Avoid 2D projection. All spatial relation tests on sphere. 

## Challenge: 
Both the points and the polygon are on the surfaces of a sphere (rather than in a cartesian space), which means that the nodes N are great circles rather than rhumb lines.

Point-in-polygon tests on a sphere are similar to point-in-polygon tests in cartesian space ([Locating a point on a spherical surface relative to a spherical polygon of arbitrary shape](http://doi.org/10.1007/BF00894449)), but appear to be more computationally expensive. PostGIS implements spherical point-in-polygon tests through [geographies](https://postgis.net/workshops/postgis-intro/geography.html).
To make it performantly feasible, we want to reduce the search space by cropping the set of points to candidates points prior to the spherical point-in-polygon tests.

# Approach:
There seems to be an opportunity to quickly retrieve candidate points through intersects test with the spherical convex hull of the polygon:

* The edges of a spherical convex hull are all great circles. 
* A great circle can be seen as a plane dividing the sphere into two hemispheres.
* A point is within the convex hull if it is in each of the hemispheres defined by each convex edge.

### Brute Force
A brute-force approach to retrieve the spherical convex hull is (as described on ([stackoverflow](https://stackoverflow.com/a/60958182):

* Consider all nodes of the polygon.
* Draw edges between each node. I.e. create a great circle that crosses each pair of nodes. Note: Draw them in both directions
* For each of the just created great circles, verify if all nodes are on the hemisphere defined by the great circle. If yes, this great circle is an edge/constraint of the convex hull. If no, discard the great circle.

### Adapted Graham Scan 
We here implement an improved approach through a spherical adaptation of the [Graham Scan](https://en.wikipedia.org/wiki/Graham_scan):

1. Find a first set of nodes that are a an edge of the convex hull. Since we know that at least one of the polygons edges is also an edge of the convex, we can look for the first convex edge in the polygon’s edges.
2. The convex hull edges have a direction; going FROM a node TO a node. 
3. We declare the first TO node as the next FROM node.
4. For this FROM node, we find its TO node. We do this by scanning all of the polygon nodes (except for the ones that already have been declared a TO nodes; but this is merely an optimization). We can discard a candidate convex edge as soon as we find a single point that is outside of its hemisphere.
5. We repeat step 3 and 4 until the TO node equals our very first FROM node.

### Further improvements through sorting:
Heuristically, the next node is probably geographically close, which probably means close in index space. The initial sorting of the nodes therefore appears relevant.

# Usage:

# Installing:

## Manual build
    
    python3 setup.py build_ext --inplace
    python3 setup.py build --force 
