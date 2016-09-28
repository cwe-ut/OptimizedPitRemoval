Optimized Pit Removal
===================

About Optimized Pit Removal
---------------------------

Extracting hydrologic features such as stream centerlines and watershed extents from a Digital Elevation Model (DEM) typically requires first hydrologically conditioning the DEM. In this process the elevations are modified in order to clearly establish flow directions. One downside of the standard ArcGIS Fill method of removing pits is that it tends to obscure meaningful elevation data for wide areas upstream of any dam-like feature. This is especially prevalent when the terrain is flat and when working with high resolution data such as LiDAR.

The Optimized Pit Removal tool uses a combination of cut and fill to remove all undesired pits while minimizing the overall changes to the landscape. The user can choose to either minimize the absolute change in landscape elevation summed across all cells, or to minimize the net change in landscape elevation (effectively balancing cut and fill). An option is also provided to exclusively use cut.

The tool also allows users to mark specific depressions to be left unmodified by setting the lowest cell to have a value of No Data. This feature can be used to establish reservoirs as well as known drainage features such as storm sewer inlets. 

The approach used by this tool is based on work by Pierre Soille (Soille, Pierre. “Optimal removal of spurious pits in grid digital elevation models.” Water Resources Research 40, W12509 (2004): 1-9.). Some changes have been made to the published algorithm. More detailed information can be found in the Documentation folder of this project.

About the tool
--------------
This tool consists of two parts. The pit removal tool itself is an executable code written in Visual Studio C++ 2010 which takes an ASCII DEM as an input and returns a hydrologically conditioned ASCII DEM as an output. This can be used via command line independently of ArcGIS. An ArcGIS 10.1 toolbox using Python contains a set of tools designed to facilitate the process of removing pits and delineating streams.

Compatability
--------------
The tool was designed for ArcMap 10.1 and Windows 7 64-bit and has not been tested in other environments. Most hydrologic conditioning workflows also require the ArcGIS Spatial Analyst Extension.

Limitations
-----------
1. The primary limitation of the current tool is that it does not support dividing a large data set into tiles for calculation. The tool treats all border cells as outlets, which may lead to unrealistic flow paths near the DEM edges. For best results, the watershed of interest should fit entirely within the DEM boundaries. In order to correct this, an algorithm will need to be developed which can remove pits in neighboring DEMs in parallel by passing information about the border cells. This is a non-trivial research task which must be solved before the Optimized Pit Removal tool can be applied to large data sets.

2. The maximum grid size is approximately 50,000,000 cells, with a practical limit of approximately 25,000,000 cells (5 square km. at 1 meter resolution). 

3. On highly detailed DEMs, the accumulated flow paths along larger streams and rivers may show an excessive number of bends, which overestimates stream lengths during flood conditions. It may be desirable to vectorize and smooth out the delineated rivers before further use. 

4. Remember that a hydrologically conditioned DEM has been modified such that it no longer matches the raw elevation data. If the conditioned DEM is intended for use other than determining flow direction, the user should carefully consider whether the adjustments to the terrain would adversely impact results.

Suggested Improvements
----------------------
1. Develop and implement an algorithm to allow pit removal of tiled or mosaicked DEMs. 
The TauDEM tools created by Dr. Tarboton support this kind of parallel processing. The algorithms used by TauDEM can serve as a starting place for developing that required by the Optimized Pit Removal tool. http://hydrology.usu.edu/taudem/taudem5/
2. Create a set of ArcGIS Toolboxes for various versions of ArcGIS.
3. Improve performance to increase the practical grid size limit.
4. Identify and correct any bugs.
5. Update tool for more recent versions of ArcGIS.
6. It has been noted that the Step Size option has difficulty if the input uses a comma rather than period for the decimal seperator. 

Related Research
----------------
LAGO Consulting and Services has developed a similar tool called Fill Sings Plus, using different principles. More information can be found here:
http://www.lago-consulting.com/fill_sinks_plus.html
