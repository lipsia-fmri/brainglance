#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 11:29:43 2019

@author: Johannes Stelzer, MPI for Biological Cybernetics, 2019
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.cm as cm

import os
import nibabel as nib
from sklearn.cluster import AffinityPropagation
import matplotlib.textpath as tp

        

class Brainglance:
    
    def __init__(self):
        self.dict_brain_areas = {}
        self.allowed_hemisphere_str = ["l", "r", "R", "L"]
        self.all_largescale_regions = []
        self.map_idx2internalindex = None
        self.S = None
        self.list_subjects = []
        
        # drawing definitions
        self.square_width = 0.8 #width of brain area rectangle (max = 1, borderless then)
        self.square_height = 0.8 #height of brain area rectangle (max = 1, borderless then)
        self.fontsize_region = 16 
        self.fontsize_area = 7
        self.edgecolor = "gray"
        self.edgewidth = 0.4 
        
        #defs spacer aka white spaces between blocks
        self.yspace_brainarealabel = 0.3
        self.xspace_region_area = 0.3
        
        self.yspace_blockedsubjects = 0.5
        self.do_mirror_lr = True #mirror the sorting of brain regions 
        self.rotregion_offset = 30 #this is the angle of rotation for the name_region text
        
        self.y_offset_cmaps = 0 #todo
        
        self.cluster_assignments = None
        self.cluster_assignment_scaling = False #scales the cluster subjects according to how many subjects are assigned to any cluster
        
        self.init_two_colormaps()
        self.do_demean = False
        
        self.dpi = 300
        
    
    def init_two_colormaps(self):
        self.set_colormaps_posneg(cm.get_cmap("Oranges"), cm.get_cmap("Blues"))
    
    def init_one_colormap(self):
        self.set_colormaps_posneg(cm.get_cmap("Oranges"))

    
    def add_atlas_definition_area(self, idx, name_area, name_largescale_region, hemisphere):
        #checks
        assert type(idx) is int, "idx needs to be an integer!"
        assert hemisphere in self.allowed_hemisphere_str, "associated_hemisphere needs to be of format {}".format(self.allowed_hemisphere_str)
        
        if idx in self.dict_brain_areas:
            print("WARNING! you have already added idx {}. Overwriting this entry...".format(idx))
        self.dict_brain_areas[idx] = [name_area, name_largescale_region, hemisphere]
        
        if name_largescale_region not in self.all_largescale_regions:
            self.all_largescale_regions.append(name_largescale_region)
            
        if self.map_idx2internalindex is None:
            self.map_idx2internalindex = idx*np.ones(1, dtype=np.int32)
        else:
            self.map_idx2internalindex = np.append(self.map_idx2internalindex, idx)
            
        self.M = self.map_idx2internalindex.size
        

    def redefine_largescale_region_sorting(self, list_sorting):
        all_largescale_regions_new = []
        for l in list_sorting:
            assert l in self.all_largescale_regions, "unknown largescale region provided: {}. Konwn largescale regions: {}".format(l, self.all_largescale_regions)
            all_largescale_regions_new.append(l)
        self.all_largescale_regions =  all_largescale_regions_new
        
        
    def add_subject(self, fp_brainmap, fp_atlas, name_subject=None, method="average"):
        assert os.path.isfile(fp_brainmap), "fp_brainmap does not exist: {}".format(fp_brainmap)
        assert os.path.isfile(fp_atlas), "fp_atlas does not exist: {}".format(fp_atlas)
        
        data_brainmap = nib.load(fp_brainmap).get_data()
        data_atlas = nib.load(fp_atlas).get_data()
        
        #atlas checks!
        assert np.issubdtype(data_atlas.dtype, np.integer), "atlas is not in integer format. Did you use the correct interpolation method? Maybe something wrong with your atlas! Make sure it is in integer representation!"
        assert data_brainmap.shape == data_atlas.shape, "brainmap and atlas size mismatch! brainmap: {}, atlas: {}".format(data_brainmap.shape, data_atlas.shape)
        
        Svec = np.zeros((1, self.M))
        for m in range(self.M):
            idx = data_atlas == self.map_idx2internalindex[m]
            data_brainmap_cut = data_brainmap[idx]
            data_brainmap_cut = data_brainmap_cut[np.isnan(data_brainmap_cut)==False]
            if data_brainmap_cut.size == 0:
                data_summary = 0
            else:
                if method == "average":
                    data_summary = np.mean(data_brainmap_cut)
                elif method == "min":
                    data_summary = np.min(data_brainmap_cut)
                elif method == "max":
                    data_summary = np.max(data_brainmap_cut)
                elif method == "std":
                    data_summary = np.std(data_brainmap_cut)
                elif method == "var":
                    data_summary = np.var(data_brainmap_cut)
                elif method == "median":
                    data_summary = np.median(data_brainmap_cut)
                    
            Svec[0, m] = data_summary
        
        if self.S is None:
            self.S = Svec
        else:
            self.S = np.append(self.S, Svec, axis=0)
            
        if name_subject is None:
            name_subject = "sub_{}".format(len(self.list_subjects))
            
        self.list_subjects.append(name_subject)
        
    def set_colormap(self, cmap):
        self.cmap = cmap
        self.use_two_colormaps = False
        
    def set_colormaps_posneg(self, cmap_pos, cmap_neg):
        self.cmap_pos = cmap_pos
        self.cmap_neg = cmap_neg
        self.use_two_colormaps = True
        
        
    def get_color(self, value):
        if self.use_two_colormaps:
            if value > 0:
                return self.cmap_pos(value)
            else:
                return self.cmap_neg(-value)      
       
    
    def get_colorbar(self, fp_figure):   
        #may only be for a subset of the data...
        idx_regionarea = []
        for r in self.all_largescale_regions:
            for lrstr in ["L","R"]:
                tmp_idx = []
                for i in range(self.S.shape[1]):
                    if self.dict_brain_areas[self.map_idx2internalindex[i]][1] == r and self.dict_brain_areas[self.map_idx2internalindex[i]][2] == lrstr:
                        tmp_idx.append(i)
                idx_regionarea.append(tmp_idx)
        all_idx_regionarea = [item for sublist in idx_regionarea for item in sublist]
        S_sub = self.S[:,all_idx_regionarea]
        maxabs = np.max(np.abs(S_sub))
        
        
        if self.use_two_colormaps:
            fig, ax = plt.subplots()
            
            fp_figure_neg = os.path.join(os.path.dirname(fp_figure), "neg_{}".format(os.path.basename(fp_figure)))
            ar=np.outer(np.ones(3),np.arange(1,0,-0.01))
            plt.imshow(ar, cmap=self.cmap_neg, extent=(0, 10, 0, 1))
            ax.set_yticks([])
            ax.set_xticks([0, 10])
            ax.set_xticklabels(["{0:2.2f}".format(np.max(self.S)), 0])
            plt.savefig(fp_figure_neg, dpi=self.dpi)
            
            fp_figure_pos = os.path.join(os.path.dirname(fp_figure), "pos_{}".format(os.path.basename(fp_figure)))
            ar=np.outer(np.ones(3),np.arange(0,1,0.01))
            plt.imshow(ar, cmap=self.cmap_pos, extent=(0, 10, 0, 1))
            ax.set_yticks([])
            ax.set_xticks([0, 10])
            ax.set_xticklabels([0, "{0:2.2f}".format(np.max(self.S))])
            plt.savefig(fp_figure_pos, dpi=self.dpi)
            
            
        else:
            fig, ax = plt.subplots()
            ar=np.outer(np.ones(3),np.arange(0,1,0.01))
            plt.imshow(ar, cmap=self.cmap_pos, extent=(0, 10, 0, 1))
            ax.set_yticks([])
            ax.set_xticks([0, 10])
            ax.set_xticklabels([0, maxabs])
            plt.savefig(fp_figure, dpi=self.dpi)

        
    def set_xspace_middle(self):
        maxwidth = 0
        for s in self.list_subjects:
            t = tp.TextPath((0,0), s, size=self.fontsize_area)
            w = t.get_extents().width
            if w > maxwidth:
                maxwidth = w
        
        self.xspace_middle = maxwidth/8
        
    def set_yspace_regiongroups(self):
        maxwidth = 0
        for key, val in self.dict_brain_areas.items():
            mx = val[0]
            t = tp.TextPath((0,0), mx, size=self.fontsize_area)
            w = t.get_extents().width
            if w > maxwidth:
                maxwidth = w
        self.yspace_regiongroups = np.cos(np.pi*self.rotregion_offset/180)*maxwidth/8
                
    
            
    def apply_clustering(self):
        #run clustering
        af = AffinityPropagation().fit(self.S) 
        cluster_centers_indices = af.cluster_centers_indices_
        labels = af.labels_   
        nmb_clusters = len(cluster_centers_indices)
        
        self.S = self.S[cluster_centers_indices,:]
        
        self.cluster_assignments = []
        for i in range(nmb_clusters):
            self.cluster_assignments.append(np.where(labels==i)[0])
            
        self.cluster_assignment_scaling = True
        self.list_subjects = ["cluster {}".format(l+1) for l in range(nmb_clusters)]
            
            
            
    def sort_subjects(self, sorting_idx=None):
        if sorting_idx is None:
            print("sorting subjects automatically (fiedler vector)")
            A = np.asmatrix(np.corrcoef(self.S))
            D = np.zeros((A.shape[0],A.shape[0]))
            L = np.asarray(D - A)
            
            try:
                evals,evec = np.linalg.eigh(L) 
                ind = np.argsort(evals) 
                evals = evals[ind] 
                evec = evec[:,ind]
                fvec = evec[:,1]
                sorting_idx = np.argsort(fvec)
            except:
                print("get_fiedlersorting problem! resorting to default sorting")
                sorting_idx = np.arange(A.shape[0])
        
        else:
            assert len(sorting_idx) == len(self.list_subjects), "sorting_idx should have the same length as the number of subjects"
            print("sorting subjects by provided sorting_idx input")
            
        
        #apply sorting
        self.list_subjects = [self.list_subjects[i] for i in sorting_idx]
        self.S = self.S[sorting_idx, :]
            
        
    def get_cluster_occupation(self, fp_figure):
        assert self.cluster_assignments is not None, "clustering was not run! run with apply_clustering."
        idx = np.arange(len(self.cluster_assignments))
        width = 0.8
        fig, ax = plt.subplots(figsize = (len(idx)/2,len(idx)/3))
        nmb_occupation = [a.size for a in self.cluster_assignments]
        ax.bar(idx, nmb_occupation, width=width)
        ax.set_xticks(idx)
        ax.set_xticklabels(self.list_subjects, rotation=90)
        plt.title("cluster occupation")
        plt.savefig(fp_figure,dpi=300)    
        

    
        

    def draw_fingerprint(self, fp_figure):
        self.set_xspace_middle()
        self.set_yspace_regiongroups()
        self.S = self.S.astype(np.float64)
        

        N = self.S.shape[0]
        M = self.S.shape[1]
        
        #convert dict into lists, first augment with hemisphere (l r)
        nmb_largescale_regions = len(self.all_largescale_regions)
        names_largescale_regions = []
        for i in range(nmb_largescale_regions):
            names_largescale_regions.append("L %s" %(self.all_largescale_regions[i]))
            names_largescale_regions.append("R %s" %(self.all_largescale_regions[i]))

        #make a list with brain region indexes (in S, not in atlas!)
        idx_regionarea = []
        names_areas = []
        for r in self.all_largescale_regions:
            for lrstr in ["L","R"]:
                tmp_idx = []
                for i in range(M):
                    if self.dict_brain_areas[self.map_idx2internalindex[i]][1] == r and self.dict_brain_areas[self.map_idx2internalindex[i]][2] == lrstr:
                        tmp_idx.append(i)
                idx_regionarea.append(tmp_idx)
                tmp_names = []
                for j in tmp_idx:
                    tmp_names.append(self.dict_brain_areas[list(self.dict_brain_areas.keys())[j]][0])
                names_areas.append(tmp_names)
        
        Mmax = max([len(m) for m in idx_regionarea]) + 1
        
        #dim stuff
        xmax = 2*Mmax+self.xspace_middle
        if self.cluster_assignment_scaling:
            nmb_subjects_percluster = np.asarray([len(l) for l in self.cluster_assignments])
            yheight_groups = np.log(nmb_subjects_percluster)+1
            yheight_groups *= 0.5
            yheight_groups_cum = np.cumsum(yheight_groups)
            yheight_groups_cum += np.arange(yheight_groups.size)*(1-self.square_height)
            ymax = yheight_groups_cum[-1]*nmb_largescale_regions+(nmb_largescale_regions)*self.yspace_regiongroups+self.y_offset_cmaps
        else:
            ymax = N*nmb_largescale_regions+(nmb_largescale_regions-1)*self.yspace_regiongroups+self.y_offset_cmaps
        
        #fig stuff
        plt.figure(figsize=(np.round(xmax/6.0),np.round(ymax/6.0)))
        plt.axis('off')
        plt.gca().set_aspect('equal')
        plt.ylim([0,ymax])
        plt.xlim((0,xmax))
        
        ax=plt.subplot(1,1,1)
        
        #normalize, so that values are between -1 and 1 (for plotting)
        #first get brain regions that are within plot scope
        all_idx_regionarea = [item for sublist in idx_regionarea for item in sublist]
        S_sub = self.S[:,all_idx_regionarea]
        
        
        if self.do_demean:
            Snorm = self.S - np.mean(self.S)
            maxabs = np.max(np.abs(Snorm))
            Snorm = Snorm/maxabs
        else:
            maxabs = np.max(np.abs(S_sub))
            Snorm = self.S/maxabs
        
        
        for g in range(nmb_largescale_regions):
            # if g<3: continue
            for lr in [0,1]:
                idx = lr+2*g
                print("Drawing largescale brainarea: {} ({}/{})".format(names_largescale_regions[idx],idx+1,len(names_largescale_regions)))
                
                #cut S matrix according to selected regiongroup
                Mi = len(idx_regionarea[idx])
                Si = Snorm[:,idx_regionarea[idx]]
                names_regiongroup = names_areas[idx]
                if lr==0: #this aligns the left hemisphere to the middle. lr=0 is the LEFT HEMISHPERE
                    Mdiff = Mmax - Mi - self.xspace_middle
                else:
                    Mdiff = 0
                
                #mirroring of LR 
                if self.do_mirror_lr and lr==0:
                    Si = Si[:,::-1]
                    names_regiongroup = names_regiongroup[::-1]
                    
                rotregion = 90 - self.rotregion_offset
                xtxtcorr = 0.5
                    
                if not self.cluster_assignment_scaling:
                    ygroupoffset = N*g + g*self.yspace_regiongroups
                    ysubj = [ymax - (s + 0.5 + self.square_height/2 + ygroupoffset) for s in range(N)]
                else:
                    ygroupoffset = yheight_groups_cum[-1]*g + g*self.yspace_regiongroups
                    ysubj = [ymax - (yheight_groups_cum[s] +ygroupoffset) for s in range(N)]
                    
                #add all squares for regiongroup
                for a in range(Mi):
                    for s in range(N):
                        cx = a + 0.5 - self.square_width/2 + Mdiff + lr*Mmax+self.xspace_middle
                        facecolor = self.get_color(Si[s,a])
                        
                        if not self.cluster_assignment_scaling:
                            boxheight = self.square_height
                        else:
                            boxheight = yheight_groups[s]
                        ax.add_patch(patches.Rectangle((cx, ysubj[s]),self.square_width, boxheight, facecolor=facecolor, edgecolor=self.edgecolor,linewidth=self.edgewidth))
                    
                    #txt individual regions 
                    y_txt = ymax + self.yspace_brainarealabel - ygroupoffset 
                    x_txt = a - self.square_width/2 + Mdiff + lr*Mmax+self.xspace_middle + xtxtcorr
                    plt.text(x_txt, y_txt, names_regiongroup[a], {'ha': 'left', 'va': 'bottom'}, fontsize=self.fontsize_area, rotation=rotregion)
                     
                    
                #txt subject subject names
                x_txt = xmax/2.0
                for s in range(N):
                    if not self.cluster_assignment_scaling:
                        y_subj_current = ysubj[s]+0.4 
                    else:
                        y_subj_current = ysubj[s]+yheight_groups[s]*0.5-0.1
                    
                    plt.text(x_txt,y_subj_current, self.list_subjects[s], {'ha': 'center', 'va': 'center'}, fontsize=self.fontsize_area)
                    
                
                #txt regiongroup names 
                if lr==0: #right hemi
                    x_txt = Mmax-Mi - self.xspace_region_area
                    bbox = {'ha': 'right', 'va': 'center'}
                else: #left hemi
                    x_txt = Mmax + Mi + self.xspace_region_area + self.xspace_middle
                    bbox = {'ha': 'left', 'va': 'center'}
                
                if not self.cluster_assignment_scaling:
                    y_txt = ymax - (np.float(N)/2  + ygroupoffset+0.1) 
                else:
                    y_txt = np.mean(ysubj) 
                plt.text(x_txt,y_txt, names_largescale_regions[idx], bbox,fontsize=self.fontsize_region)
        
        
        
        plt.savefig(fp_figure,dpi=self.dpi)
        
    

   