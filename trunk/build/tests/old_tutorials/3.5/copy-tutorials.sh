
#mkdir -p toast ; svn add toast
mkdir -p GREIT-evaluation ; svn add GREIT-evaluation
mkdir -p GREIT-evaluation/EIT08-GREIT-pres ; svn add GREIT-evaluation/EIT08-GREIT-pres
mkdir -p GREIT-evaluation/simulation_test_imgs ; svn add GREIT-evaluation/simulation_test_imgs
mkdir -p absolute_reconst ; svn add absolute_reconst
mkdir -p data_structures ; svn add data_structures
mkdir -p version2 ; svn add version2
mkdir -p GREIT ; svn add GREIT
mkdir -p vtk_visualization ; svn add vtk_visualization
mkdir -p EIDORS_basics ; svn add EIDORS_basics
mkdir -p lung_EIT ; svn add lung_EIT
mkdir -p cheating_EIDORS ; svn add cheating_EIDORS
mkdir -p dual_model ; svn add dual_model
mkdir -p GREIT-algorithm ; svn add GREIT-algorithm
mkdir -p distmesh ; svn add distmesh
mkdir -p adv_image_reconst ; svn add adv_image_reconst
mkdir -p adv_image_reconst/tv_hp_imgs ; svn add adv_image_reconst/tv_hp_imgs
mkdir -p netgen ; svn add netgen
mkdir -p netgen/extrusion ; svn add netgen/extrusion
mkdir -p geophysics ; svn add geophysics
mkdir -p strange_effects ; svn add strange_effects
mkdir -p other_models ; svn add other_models

#svn cp ../../../../../tags/3.5/htdocs/tutorial/toast/*.m                                    toast/
svn cp ../../../../../tags/3.5/htdocs/tutorial/GREIT-evaluation/*.m                         GREIT-evaluation/
svn cp ../../../../../tags/3.5/htdocs/tutorial/GREIT-evaluation/EIT08-GREIT-pres/*.m        GREIT-evaluation/EIT08-GREIT-pres/
svn cp ../../../../../tags/3.5/htdocs/tutorial/GREIT-evaluation/simulation_test_imgs/*.m    GREIT-evaluation/simulation_test_imgs/
svn cp ../../../../../tags/3.5/htdocs/tutorial/absolute_reconst/*.m                         absolute_reconst/
svn cp ../../../../../tags/3.5/htdocs/tutorial/data_structures/*.m                          data_structures/
svn cp ../../../../../tags/3.5/htdocs/tutorial/version2/*.m                                 version2/
svn cp ../../../../../tags/3.5/htdocs/tutorial/GREIT/*.m                                    GREIT/
svn cp ../../../../../tags/3.5/htdocs/tutorial/vtk_visualization/*.m                        vtk_visualization/
svn cp ../../../../../tags/3.5/htdocs/tutorial/EIDORS_basics/*.m                            EIDORS_basics/
svn cp ../../../../../tags/3.5/htdocs/tutorial/lung_EIT/*.m                                 lung_EIT/
svn cp ../../../../../tags/3.5/htdocs/tutorial/cheating_EIDORS/*.m                          cheating_EIDORS/
svn cp ../../../../../tags/3.5/htdocs/tutorial/dual_model/*.m                               dual_model/
svn cp ../../../../../tags/3.5/htdocs/tutorial/GREIT-algorithm/*.m                          GREIT-algorithm/
svn cp ../../../../../tags/3.5/htdocs/tutorial/distmesh/*.m                                 distmesh/
svn cp ../../../../../tags/3.5/htdocs/tutorial/adv_image_reconst/*.m                        adv_image_reconst/
svn cp ../../../../../tags/3.5/htdocs/tutorial/adv_image_reconst/tv_hp_imgs/*.m             adv_image_reconst/tv_hp_imgs/
svn cp ../../../../../tags/3.5/htdocs/tutorial/netgen/*.m                                   netgen/
svn cp ../../../../../tags/3.5/htdocs/tutorial/netgen/extrusion/*.m                         netgen/extrusion/
svn cp ../../../../../tags/3.5/htdocs/tutorial/geophysics/*.m                               geophysics/
svn cp ../../../../../tags/3.5/htdocs/tutorial/strange_effects/*.m                          strange_effects/
svn cp ../../../../../tags/3.5/htdocs/tutorial/other_models/*.m                             other_models/
