<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xsl" href="OpenGeoSysGLI.xsl"?>
<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogs="http://www.opengeosys.org">
    <name>beacon_1c_geometry</name>
    <points>
        <point id="0" x="0.0" y="0.0" z="0.0" name="origin"/>
        <point id="1" x="0.0" y="0.0485" z="0.0"/>
        <point id="2" x="0.0" y="0.1" z="0.0"/>
        <point id="3" x="0.05" y="0.0" z="0.0"/>
        <point id="4" x="0.05" y="0.0485" z="0.0"/>
        <point id="5" x="0.05" y="0.1" z="0.0"/>
    </points>
    <polylines>
        <polyline id="0" name="left"><pnt>0</pnt><pnt>1</pnt><pnt>2</pnt></polyline>
        <polyline id="1" name="right"><pnt>3</pnt><pnt>4</pnt><pnt>5</pnt></polyline>
        <polyline id="2" name="bottom"><pnt>0</pnt><pnt>3</pnt></polyline>
        <polyline id="3" name="top"><pnt>2</pnt><pnt>5</pnt></polyline>
    </polylines>
</OpenGeoSysGLI>
