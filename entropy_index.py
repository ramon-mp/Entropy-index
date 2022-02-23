# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 11:21:55 2022
@author: Molinero-Parejo, R.
"""

from os import chdir
from geopandas import read_file
from time import time
import numpy as np

def entropy_index(parcel, block):
    """
    Function that will calculate the entropy index for each block
    """  
    
    # generate new gdf
    gdf_p = parcel[{"REFCAT", "USO", "AREA", "geometry"}]
    gdf_b = block[{"MASA", "AREA", "geometry"}]
    
    # generate new column MASA to join "one to many"
    gdf_p["MASA"] = gdf_p["REFCAT"].apply(lambda x: x[0:5])
    
    # generate new df contain area of each block (MASA)
    gdf_area_m = gdf_p[["MASA", "AREA"]].groupby("MASA").sum("AREA")
    
    # merge parcels with blocks area
    gdf_merge = gdf_p.merge(gdf_area_m, how="inner", on = "MASA", suffixes=("_P", "_M"))
    
    # list of active uses
    list_uses = (gdf_p["USO"].unique()).tolist()
    list_uses.remove("SIN_EDIF")
    
    # create one field for each use
    for use in list_uses:
        gdf_merge[use] = 0
    
    # calculate proportion of each parcel inside a block
    for index, row in gdf_merge.iterrows():
        for use in list_uses:
            if use == row["USO"]:
                gdf_merge.loc[index, use] = row["AREA_P"] / row["AREA_M"]
            else:
                gdf_merge.loc[index, use] = 0
    
    # group parcels of the same use inside a block
    list_merge_fields = ["MASA", "USO", "AREA_P"] + list_uses
    list_key = ["AREA_P"] + list_uses
    list_value = ["sum"] * len(list_key)
    dict_uses = dict(zip(list_key, list_value))
    gdf_prop = gdf_merge[list_merge_fields].groupby(["MASA", "USO"]).agg(dict_uses).reset_index()

    # calculate LN for the metric formula -[P * LN(P)] - NUMERADOR
    for use in list_uses:
        gdf_prop["LN_" + use] = gdf_prop[use].apply(lambda x: np.log(x) if x != 0 else 0)
        gdf_prop["X_" + use] = gdf_prop[use] * gdf_prop["LN_" + use]
    list_x_uses = []
    for use in list_uses:
        list_x_uses.append("X_" + use)    
    gdf_prop["NUMERADOR"] = gdf_prop[list_x_uses].sum(axis=1)
    
    # delete SIN_EDIF from uses
    gdf_prop.loc[gdf_prop.USO == "SIN_EDIF", "USO"] =  float("NaN")

    # group rows by block ID (MASA)
    gdf_entropy = gdf_prop[["MASA", "NUMERADOR", "USO"]].groupby("MASA").agg({"NUMERADOR": "sum", "USO": "nunique"}).reset_index()
    
    # calculate LN for the metric formula LN(k) - DENOMINADOR
    gdf_entropy["DENOMINADOR"] = gdf_entropy["USO"].apply(lambda x: np.log(x))
    
    # calculate formula -[P * LN(P)] / LN(k)
    gdf_entropy["ENTROPY"] = -(gdf_entropy["NUMERADOR"]) / gdf_entropy["DENOMINADOR"]
    
    # fill Nan and Inf with 0
    gdf_entropy = gdf_entropy.fillna(0)
    gdf_entropy = gdf_entropy.replace([np.inf, -np.inf], 0)
    
    # join parcel df to block df and export (*.shp)
    gdf_b_index = gdf_b.merge(gdf_entropy, how="inner", on = "MASA", suffixes=("_0", "_1"))
    gdf_b_index.to_file("block_mix.shp")
    
#-----------------------------------------------------------------------------#

def main():
    """
    Main function that will set up the data and organize it so that the 
    processes are automatized
    """
    
    try:
        # open PARCEL file as GeoDataFrame
        shp_parcel = "clas_i_urb.shp"
        gdf_parcel = read_file(shp_parcel)
        
        # open MASA file as GeoDataFrame
        shp_block = "MASA.shp"
        gdf_block = read_file(shp_block)
    
        print("Opening file: " + shp_parcel + "\n")
        print("Opening file: " + shp_block + "\n")
    
        # run the function
        entropy_index(gdf_parcel, gdf_block)

        print("Saving file")

    except ValueError:
        print("\n---> ERROR: " + ValueError)

#-----------------------------------------------------------------------------#

# save initial time
t_start = time()  

# establish the working directory
wd = "D:\\pruebas"
chdir(wd)

print("Executing the script." + "\n")
print("Setting working directory in: " + wd + "\n")

# executing main function
if __name__ == "__main__":
    main()

# save the final time
t_finish = time()

# show the execution time elapsed
t_process = (t_finish - t_start) / 60
print("Process time: " + str(round(t_process, 2)) + "minutes")