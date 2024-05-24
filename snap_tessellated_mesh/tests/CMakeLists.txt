include(CTest)
add_test(NAME RunCurvedWall COMMAND 
         snap_to_tessellated_mesh.exe
           ${CMAKE_CURRENT_SOURCE_DIR}/curved_wall_input.vtk
           ${CMAKE_CURRENT_SOURCE_DIR}/curved_wall_input.igs 
           test_output.vtk)
add_test(NAME CompareCurvedWall COMMAND 
         ${CMAKE_COMMAND} -E compare_files
         test_output.vtk ${CMAKE_CURRENT_SOURCE_DIR}/curved_wall_output.vtk)

add_test(NAME RunCurvedWallNoZFaces COMMAND
         snap_to_tessellated_mesh.exe
           ${CMAKE_CURRENT_SOURCE_DIR}/curved_wall_input.vtk
           ${CMAKE_CURRENT_SOURCE_DIR}/curved_wall_input.igs
           test_output_z.vtk
           1)
add_test(NAME CompareCurvedWallNoZFaces COMMAND
         ${CMAKE_COMMAND} -E compare_files
         test_output_z.vtk ${CMAKE_CURRENT_SOURCE_DIR}/curved_wall_output_z.vtk)
