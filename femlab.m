% FEMLAB
%
% This is a set of Matlab routines do perform some basic finite element
% routines.  Currently this consistes of teh following
%
% Mesh Generation Routines:
%   readelements.m
%   readelemstruct.m
%   readnodes.m
%   readnodeset.m
%   renumber.m
%
%   nodearray1d
%   nodearray2d
%   nodearray3d
%
%   genconn2d
%   genconn3d
%
%   plot_fefield
%   plot_mesh
%
% Shape Function Routines:
%   shape.m
%   shape_hexa8.m
%   shape_line2.m
%   shape_line3.m
%   shape_quad4.m
%   shape_tetra4.m
%   shape_tria3.m
%
%   dshape.m
%   dshape_hexa8.m
%   dshape_line2.m
%   dshape_line3.m
%   dshape_quad4.m
%   dshape_tetra4.m
%   dshape_tria3.m
%
% Element B-matricies:
%   bmat_beam2d.m
%   bmat_beam3d.m
%   bmat_hexa20.m
%   bmat_hexa8.m
%   bmat_quad4.m
%   bmat_quad8.m
%   bmat_tetra10.m
%   bmat_tetra4.m
%   bmat_tria3.m
%   bmat_tria6.m
%   bmat_truss2d.m
%   bmat_truss3d
%   bmat_truss3d.m
%
% Material defs:
%   cmat_mat1.m
%
% Basic Element Formulation:
%   element_jacobian.m
%   fill_bmatrix.m
%   get_scatter.m
%   ShellDirectionCosines.m
%   grad_shapefunct.m
%   element_porder.m
%   etypestr.m
%   element_quadrature.m
%   quadrature_gaussian.m
%   quadrature_simplex.m
%   quadrature_compound.m
%
% K-Matrix computations (linear elasticity)
%   kmat_beam2d.m
%   kmat_beam3d.m
%   kmat_hexa20.m
%   kmat_hexa8.m
%   kmat_mindlinshell.m
%   kmat_quad4.m
%   kmat_quad8.m
%   kmat_tetra10.m
%   kmat_tetra4.m
%   kmat_tria3.m
%   kmat_tria6.m
%   kmat_truss2d.m
%   kmat_truss3d.m
% 
% Load calculations:
%   uniform_trac2d.m
%
% Basic Postprocessing:
%   nodal_avg.m
%   principal_val.m
%
% Ensight/Paraview output:
%   ensight_case.m
%   ensight_fegeometry.m
%   ensight_field.m
%
% Solving:
%   fesolve.m
%   fesolve2.m
%   rotkdof.m
%
% Examples: 
%   fea1d.m
%   fea2d.m
%   fea2d_tria3.m
%   fea2d_truss.m
%   heat_conduction.m
%   test1.m
%

% Examples
% LevelSet
% ShellDirectionCosines.m
% XFEM
% bmat_axisymmetric.m
% bmat_beam2d.m
% bmat_beam3d.m
% bmat_hexa20.m
% bmat_hexa8.m
% bmat_quad4.m
% bmat_quad8.m
% bmat_tetra10.m
% bmat_tetra4.m
% bmat_tria3.m
% bmat_tria6.m
% bmat_truss2d.m
% bmat_truss3d.m
% cmat_mat1.m
% dshape.m
% dshape_hexa8.m
% dshape_line2.m
% dshape_line3.m
% dshape_quad4.m
% dshape_tetra4.m
% dshape_tria3.m
% dshape_tria6.m
% element_jacobian.m
% element_jacobian2.m
% element_porder.m
% element_quadrature.m
% ensight_blockgeometry.m
% ensight_case.m
% ensight_fegeometry.m
% ensight_field.m
% etypestr.m
% extrude_nodes.m
% femlab.m
% fesolve.m
% fesolve2.m
% fill_bmatrix.m
% genconn2d.m
% genconn3d.m
% get_scatter.m
% grad_shapefunct.m
% grad_shapefunct2.m
% kmat_BT_shell.m
% kmat_beam2d.m
% kmat_beam3d.m
% kmat_hexa20.m
% kmat_hexa8.m
% kmat_mindlin_quad4.m
% kmat_quad4.m
% kmat_quad8.m
% kmat_tetra10.m
% kmat_tetra4.m
% kmat_tria3.m
% kmat_tria6.m
% kmat_truss2d.m
% kmat_truss3d.m
% meshing
% mises_val.m
% nodal_avg.m
% nodearray1d.m
% nodearray2d.m
% nodearray3d.m
% plot_fefield.m
% plot_mesh.m
% principal_val.m
% quadrature_compound.m
% quadrature_gaussian.m
% quadrature_simplex.m
% read_dyna3d.m
% readelements.m
% readelemstruct.m
% readnodes.m
% readnodeset.m
% renumber.m
% rotkdof.m
% shape.m
% shape_axisymmetric.m
% shape_cubic_hermite.m
% shape_hexa8.m
% shape_line2.m
% shape_line3.m
% shape_quad4.m
% shape_tetra4.m
% shape_tria3.m
% shape_tria6.m
% uniform_trac2d.m