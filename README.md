BRAINGLANCE
================================
Brainglance is a python-based software package to display single subject information for MRI.

![brainglance_example](https://github.com/lipsia-fmri/brainglance/blob/master/example.png)






INSTALLATION
================================

To install brainglance, make sure your python distribution satisfies all requirements listed.
You can install them by



        pip install -r requirements.txt


BASIC USAGE
================================

PREREQUISITES
-----------------

You will need the following inputs:

- a 3D brain map of each subject in native subject space, containing some kind of data you want to display
- a coregistered atlas in native subject space (e.g. Brainnetome, AAL, ...). Make sure you have used nearest neighbor interpolation here, as only integer values are accepted.
- a description of the atlas, where the label is associated with the name of the brain region, the name of the large scale region and the hemisphere. For instance label=1, name_region='SFG_7_1', name_largescale_region='FRO', hemisphere='L'. This type of information usually is supplied with every atlas.


RUNNING BRAINGLANCE
--------------------------

The first step is to instantiate the brainglance object:


        bg = BrainGlance()




You will need to fill the description of the atlas. The "label" is corresponding to the integer value in the atlas. You can add all information manually, e.g. manually:

        bg.add_atlas_definition_area(label=1, name_area="area1", name_largescale_region="FRO", hemisphere="L")
        bg.add_atlas_definition_area(label=2, name_area="area2", name_largescale_region="FRO", hemisphere="L")
        bg.add_atlas_definition_area(label=3, name_area="area1", name_largescale_region="FRO", hemisphere="R")
        bg.add_atlas_definition_area(label=4, name_area="area2", name_largescale_region="FRO", hemisphere="R")


In case of the Brainnetome atlas, the following code achieves this:


        M = 246
        dp_atlas='/media/3tbd/studies/_data'
        fn_atlas_descr = os.path.join(dp_atlas,'BNA_brainatlas_areas.txt')
        atlas_descr = np.loadtxt(fn_atlas_descr,dtype=str,comments='#',delimiter='\t')
        atlas_descr = atlas_descr[1:M+1,:]


Next, we will read in this list into the brainglance object.



        for i in range(M):
            label = np.int(atlas_descr[i][0])
            name_area = atlas_descr[i][1]
            name_largescale_region = atlas_descr[i][5].upper()
            hemisphere = atlas_descr[i][2].upper()
            bg.add_atlas_definition_area(label, name_area, name_largescale_region, hemisphere)


Now you need to add your subjects. For each subject, you need to provide the 3D brain map with values you want to show (fp_brainmap) and a coregistered atlas in native subject space (fp_atlas). Furthermore, you can supply the name of the subject as third argument, however, this is optional.

Add a single subject via


          bg.add_subject(fp_brainmap, fp_atlas, "subject-01")

You can add more subjects with the same call:

          bg.add_subject(fp_brainmap2, fp_atlas2, "subject-02")
          bg.add_subject(fp_brainmap3, fp_atlas3, "subject-03")


Now you can generate a brainglance plot, supplying the figure output file fp_figure.


          bg.draw_fingerprint(fp_figure)



EXTRA OPTIONS
=========================

Change order or selection of brain regions
-----------------------------------------------

If you want to change the order of the largescale regions or only plot a subset of largescale regions, you can supply a list with the names (and order). In the following, we want to restrict the plotting to temporal (TEM) and frontal (FRO) regions


          bg.redefine_largescale_region_sorting(["TEM", "FRO"])
          bg.draw_fingerprint(fp_figure)


Change sorting of subjects
-----------------------------------------------

There are two ways of sorting the subjects: either you supply a list with the sorting, e.g. sorting_idx = [1, 0, 3, ...], or alternatively, use automatic sorting according to the Fielder vector.


          bg.sort_subjects(sorting_idx) #manual sorting
          bg.sort_subjects() #automatic fiedler sorting



Set colormaps
------------------
You can set either one or two colormaps. In the latter case, one colormap will be used for positive values and one for negative ones. You will need to pass the colormap defined in matplotlib



         import matplotlib.cm as cm
         cmap = cm.get_cmap("Oranges")
         bg.set_colormap(cmap) #one colormap
         bg.set_colormaps_posneg(cmap_pos, cmap_neg) #two colormaps



Set demeaning
------------------

If you want to demean your data for visualization, set bg.do_demean=True. This can be useful for data that is strictly positive, however you want to visualize it with two colormaps.


Get colorbars
--------------------

To print the colorbars, simply supply a filename fp_figure and run:


         bg.get_colorbar(fp_figure)




Apply clustering
-----------------------------------------------

For clustering the subject into prototypes and displaying these weighted by their respective occupations, simply call:


          bg.apply_clustering()
          bg.draw_fingerprint(fp_figure)


Note: this changes the underlying data structures. If you want to switch back to displaying single subjects, you will need to re-initialize the bg object from scratch.


Print occupation of clusters
-----------------------------------------------

To create a histogram of the cluster occupatio and save it with the filename fp_figure, simply run:


          bg.get_cluster_occupation(fp_figure)
