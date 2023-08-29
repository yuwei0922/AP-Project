# Photogrammetry
A C++ project that implements the five core algorithms of photogrammetry (analytical part), including spatial resection, spatial intersection, interior orientation, relative orientation, and absolute orientation.

1. Spatial Resection:
Input: Corresponding pairs of image points and ground coordinates
Output: Six exterior orientation elements, rotation matrix R

2. Forward Intersection:
Methods of point projection coefficients and rigorous solution of collinearity equations
Input: Image point coordinates, camera parameters, exterior orientation elements
Output: Ground coordinates computed using both methods

3. Interior Orientation:
Input: Coordinates of fiducial marks, corresponding image point coordinates
Output: Interior orientation parameters

4. Relative Orientation:
Input: Focal length, image point coordinates
Output: Relative orientation elements, residuals for each image point, model point coordinates

5. Absolute Orientation:
Input: Corresponding pairs of model points and ground coordinates
Output: Seven-parameter absolute orientation solution, ground photogrammetric coordinates for the unknown points
