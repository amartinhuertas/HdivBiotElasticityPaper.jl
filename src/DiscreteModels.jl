function generate_model_unit_square_biot_elasticity(nk;
                                                      simplexify_model=true,
                                                      bcs_type=:mixed)
    domain    = (0.0,1.0,0.0,1.0)
    n         = 2^nk
    partition = (n,2*n)
    m=generate_rectangular_biot_elasticity_mesh(:t,
                                                 n,
                                                 [:b],
                                                 domain,
                                                 partition,
                                                 bcs_type=bcs_type)
    if (simplexify_model)
      m=simplexify(m)
    end
    m
end

function set_boundary_vertices_entity_to_facets_entity!(model)
  labeling = get_face_labeling(model)
  topology = Gridap.Geometry.get_grid_topology(model)
  D = num_cell_dims(topology)
  d = D-1
  face_to_cells    = Gridap.Geometry.get_faces(topology,d,D)
  face_to_vertices = Gridap.Geometry.get_faces(topology,d,0)
  face_to_cells_ptrs = face_to_cells.ptrs
  nfaces=length(face_to_cells_ptrs)-1
  for face=1:nfaces
    ncells_around = face_to_cells_ptrs[face+1] - face_to_cells_ptrs[face]
    if ncells_around == 1
      face_entity = labeling.d_to_dface_to_entity[d+1][face]
      s=face_to_vertices.ptrs[face]
      e=face_to_vertices.ptrs[face+1]-1
      for j=s:e
        vertex=face_to_vertices.data[j]
        labeling.d_to_dface_to_entity[1][vertex]=face_entity
      end
    end
  end
end

function set_vessel_interstitium_cell_entities!(face_labeling,
                                                NX,
                                                NY,
                                                vessel_ncells_in_thin_dim,
                                                vessel_entity,
                                                inter_entity,
                                                vessel_location) 
    face_labeling.d_to_dface_to_entity[3] .= inter_entity
    lid=LinearIndices((NX,NY))
    for vl in vessel_location
        if (vl==:b)
            s=1
            e=vessel_ncells_in_thin_dim*NX
            face_labeling.d_to_dface_to_entity[3][s:e] .= vessel_entity
        elseif (vl==:t)
            e=length(face_labeling.d_to_dface_to_entity[3])
            s=NX*NY-vessel_ncells_in_thin_dim*NX+1
            face_labeling.d_to_dface_to_entity[3][s:e] .= vessel_entity
        elseif (vl==:l)
            for j=1:NY
                for i=1:vessel_ncells_in_thin_dim
                    face_labeling.d_to_dface_to_entity[3][lid[i,j]] = vessel_entity
                end
            end
        elseif (vl==:r)
            for j=1:NY
                for i=NX-vessel_ncells_in_thin_dim+1:NX
                    face_labeling.d_to_dface_to_entity[3][lid[i,j]] = vessel_entity
                end
            end
        end
    end
end

function set_face_entities!(face_labeling,
                            grid_topology,
                            NX,
                            NY,
                            vessel_entity,
                            inter_entity,
                            vessel_face_entity,
                            inter_face_entity,
                            boundary_location)
  @assert boundary_location in (:t,:b,:l,:r)
  cell_faces=Gridap.Geometry.get_faces(grid_topology,2,1)
  cell_vertices=Gridap.Geometry.get_faces(grid_topology,2,0)
  lid=LinearIndices((NX,NY))
  if (boundary_location==:b)
    verts=[1,2]
    lfid=1
    j=1
    for i=1:NX
      cell=lid[i,j]
      if face_labeling.d_to_dface_to_entity[3][cell] == inter_entity
        face_labeling.d_to_dface_to_entity[2][cell_faces[cell][lfid]] = inter_face_entity
        face_entity=inter_face_entity
      elseif face_labeling.d_to_dface_to_entity[3][cell] == vessel_entity
        face_labeling.d_to_dface_to_entity[2][cell_faces[cell][lfid]] = vessel_face_entity
        face_entity=vessel_face_entity
      else 
        @assert false  
      end
      for v in verts 
        face_labeling.d_to_dface_to_entity[1][cell_vertices[cell][v]] = face_entity
      end
    end
  elseif (boundary_location==:t)
    verts=[3,4]
    lfid=2
    j=NY
    for i=1:NX 
        cell=lid[i,j]
        if face_labeling.d_to_dface_to_entity[3][cell] == inter_entity
            face_labeling.d_to_dface_to_entity[2][cell_faces[cell][lfid]] = inter_face_entity
            face_entity=inter_face_entity
        elseif face_labeling.d_to_dface_to_entity[3][cell] == vessel_entity
            face_labeling.d_to_dface_to_entity[2][cell_faces[cell][lfid]] = vessel_face_entity
            face_entity=vessel_face_entity
        else 
            @assert false  
        end
        for v in verts 
          face_labeling.d_to_dface_to_entity[1][cell_vertices[cell][v]] = face_entity
        end
    end
  elseif (boundary_location==:l)
    verts=[1,3]
    lfid=3
    i=1
    for j=1:NY
        cell=lid[i,j]
        if face_labeling.d_to_dface_to_entity[3][cell] == inter_entity
            face_labeling.d_to_dface_to_entity[2][cell_faces[cell][lfid]] = inter_face_entity
            face_entity=inter_face_entity
        elseif face_labeling.d_to_dface_to_entity[3][cell] == vessel_entity
            face_labeling.d_to_dface_to_entity[2][cell_faces[cell][lfid]] = vessel_face_entity
            face_entity=vessel_face_entity
        else 
            @assert false  
        end
        for v in verts 
          face_labeling.d_to_dface_to_entity[1][cell_vertices[cell][v]] = face_entity
        end
    end
  elseif (boundary_location==:r)
    verts=[2,4]
    lfid=4
    i=NX
    for j=1:NY
        cell=lid[i,j]
        if face_labeling.d_to_dface_to_entity[3][cell] == inter_entity
            face_labeling.d_to_dface_to_entity[2][cell_faces[cell][lfid]] = inter_face_entity
            face_entity=inter_face_entity
        elseif face_labeling.d_to_dface_to_entity[3][cell] == vessel_entity
            face_labeling.d_to_dface_to_entity[2][cell_faces[cell][lfid]] = vessel_face_entity
            face_entity=vessel_face_entity
        else 
            @assert false  
        end
        for v in verts 
            face_labeling.d_to_dface_to_entity[1][cell_vertices[cell][v]] = face_entity
        end
    end
  end
end

function generate_rectangular_biot_elasticity_mesh(block,
                                                   ncells_elasticity_region_thin_dim,
                                                   E_location::Vector{Symbol},
                                                   domain,
                                                   partition;
                                                   bcs_type=:mixed,
                                                   Etag="elast",
                                                   Ptag="poroelast",
                                                   wallEDtag="wallED",
                                                   wallPDtag="wallPD",
                                                   wallENtag="wallEN",
                                                   wallPNtag="wallPN")
  @assert bcs_type in (:mixed,:full_dirichlet) 

  model=CartesianDiscreteModel(domain,partition)
  NX,NY=partition
  face_labeling=model.face_labeling
  max_cell_entity=maximum(face_labeling.d_to_dface_to_entity[3])
  E_entity=max_cell_entity+1
  P_entity=max_cell_entity+2
  Gridap.add_tag!(face_labeling,Etag,[E_entity])
  Gridap.add_tag!(face_labeling,Ptag,[P_entity])
  set_vessel_interstitium_cell_entities!(face_labeling,
                                          NX,
                                          NY,
                                          ncells_elasticity_region_thin_dim,
                                          E_entity,
                                          P_entity,
                                          E_location)

  grid_topology = Gridap.Geometry.get_grid_topology(model)

  EN_face_entity = max_cell_entity+5
  PN_face_entity = max_cell_entity+6
  if (block == :t)
        for i in (:l,:r)
            set_face_entities!(face_labeling,
                            grid_topology,
                            NX,
                            NY,
                            E_entity,
                            P_entity,
                            EN_face_entity,
                            PN_face_entity,
                            i)
        end 
  elseif (block == :bl || block == :br)
        for i in (:l,:r,:t)
            set_face_entities!(face_labeling,
                            grid_topology,
                            NX,
                            NY,
                            E_entity,
                            P_entity,
                            EN_face_entity,
                            PN_face_entity,
                            i)
        end 
  end

  ED_face_entity = max_cell_entity+3
  PD_face_entity  = max_cell_entity+4

  if (bcs_type==:mixed)
    Gridap.add_tag!(face_labeling,wallENtag,[EN_face_entity])
    Gridap.add_tag!(face_labeling,wallPNtag,[PN_face_entity])
    Gridap.add_tag!(face_labeling,wallEDtag,[ED_face_entity])
    Gridap.add_tag!(face_labeling,wallPDtag,[PD_face_entity])
  else
    Gridap.add_tag!(face_labeling,wallENtag,Int[])
    Gridap.add_tag!(face_labeling,wallPNtag,Int[])
    Gridap.add_tag!(face_labeling,wallEDtag,[ED_face_entity,EN_face_entity])
    Gridap.add_tag!(face_labeling,wallPDtag,[PD_face_entity,PN_face_entity])
  end 

  if (block == :t)
    set_face_entities!(face_labeling,
                        grid_topology,
                        NX,
                        NY,
                        E_entity,
                        P_entity,
                        ED_face_entity,
                        PD_face_entity,
                        :t)
    set_face_entities!(face_labeling,
    grid_topology,
    NX,
    NY,
    E_entity,
    P_entity,
    ED_face_entity,
    PD_face_entity,
    :b)
  elseif (block == :bl || block == :br)
    set_face_entities!(face_labeling,
                        grid_topology,
                        NX,
                        NY,
                        E_entity,
                        P_entity,
                        ED_face_entity,
                        PD_face_entity,
                        :b)
  end 
  model
end 


function generate_vessel_interstitium_model_block(block,
                                                  vessel_ncells_in_thin_dim,
                                                  vessel_location::Vector{Symbol})
    @assert block in (:t,:bl,:br)
    @assert length(vessel_location)>=1 && length(vessel_location)<=2
    @assert all([x in (:t,:b,:l,:r) for x in vessel_location])
    vessel_thickness = 0.01
    h_vessel = vessel_thickness/vessel_ncells_in_thin_dim
    if block==:t 
        ymin = 0.17
        ymax = 0.25
        xmin = 0.0
        xmax = 0.25
      elseif block==:bl
        ymin = 0.0
        ymax = 0.13
        xmin = 0.0 
        xmax = 0.105
      elseif block==:br
        xmin = 0.145
        xmax = 0.25 
        ymin = 0.0 
        ymax = 0.13  
      end 
      NY = round(Int,(ymax - ymin)/h_vessel)
      NX = round(Int, (xmax - xmin)/h_vessel)
      @assert abs(vessel_ncells_in_thin_dim*(ymax - ymin)/NY) ≈ vessel_thickness
      domain=(xmin,xmax,ymin,ymax)
      partition=(NX,NY)
      generate_rectangular_biot_elasticity_mesh(block,
                                                vessel_ncells_in_thin_dim,
                                                vessel_location,
                                                domain,
                                                partition,
                                                bcs_type=:mixed,
                                                Etag="vessel",
                                                Ptag="interstitium",
                                                wallEDtag="clamped_vessel",
                                                wallPDtag="clamped_interstitium",
                                                wallENtag="free_vessel",
                                                wallPNtag="free_interstitium")
end