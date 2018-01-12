"""
Helper functions to manipulate xml model files
"""

# --- Imports ----- #
import logging
from lxml.etree import Element, SubElement, parse
# ----------------- #

def change_file_function(xmlfile, filefunction, srcname):
    """
    Read in an xml model file and change the file function for a particular source.
    Overwrites the existing xmlfile.

    Parameters
    ----------
    xmlfile:	string,
		path to xml model file
    filefunction:	string, 
		path to new file function / map cube fits file
    srcname:	string,
		name of the source for which file function should be changed
    """
    f = open(xmlfile)
    tree = parse(f)
    f.close()
    sourceLib = tree.getroot()
    for src in sourceLib:
	if not src.get('name') == srcname: continue
	for model in src:
	    if not model.tag == 'spatialModel': continue
	    if not model.get('type') == 'MapCubeFunction':
		raise ValueError("Spatial model must be of type MapCubeFunction, is {0:s}".format(
		    model.get('type')))
	    else:
		model.set("file", filefunction)
		logging.info("Set file for src {0:s} to {1:s}".format(src.get('name'),
		    model.get('file')))
    tree.write(xmlfile)
    logging.info("Overwritten XML file {0:s}".format(xmlfile))
    return 

def add_src_file_function(xmlfile, filefunction, srcname, ra = 0., dec = 0.):
    """
    Read in an xml model file and add a new point source model with a filefunction.
    Overwrites the existing xmlfile.

    Parameters
    ----------
    xmlfile:	string,
		path to xml model file
    filefunction:	string, 
		path to new file function
    srcname:	string,
		name of the new point source 
    ra: float
    		Right ascension of new source in degrees
    dec: float
    		Right ascension of new source in degrees
    """
    f = open(xmlfile)
    tree = parse(f)
    f.close()
    sourceLib = tree.getroot()
    for src in sourceLib:
	for model in src:
	    if not model.tag == 'spatialModel': continue
	    if model.get('type') == 'SkyDirFunction':
		new_src = src # copy an exisiting point source
    new_src.set('name', srcname)
    for model in new_src:
	if model.tag == 'spatialModel': 
	    for p in model:
		if p.get('name') == 'RA':
		    p.set('value',str(ra))
		elif p.get('name') == 'DEC':
		    p.set('value',str(dec))
	if model.tag == 'spectrum': 
	    model.set('type','FileFunction')
	    model.set('file',filefunction)
	    if len(model) > 1:
		for t in model[1:]:
		    model.remove(t)
	    params = dict(free="1", max="1e5", min="1e-5", name="Normalization", scale="1", value="1.0")
	    for k,v in params.items():
		model[0].set(k, v)

    logging.info("Added new source {0:s} with FileFunction {1:s} at coordinates ({2:.2f},{3:.2f})".format(
    		srcname, filefunction, ra, dec))
    tree.write(xmlfile)
    logging.info("Overwritten XML file {0:s}".format(xmlfile))
    return

