PPRINT/ --- files_x\StartEnd\header.txt ---
$$ HYPERMILL
$$*
PPRINT/ --- files_x\StartEnd\file_start.txt ---
$$ -> MFGNO
PARTNO / Mold Repair 1 1 Main Spindle Machining
MACHIN / MILL, 01
$$ -> CUTCOM_GEOMETRY_TYPE / 
UNITS / MM
CALSUB/START_PROG
PARTNO/1 1 Main Spindle Machining
PPRINT/ --- files_x\ToolChange\Toollist_Begin.txt ---
PPRINT/ --- files_x\ToolChange\Toollist_End.txt ---
PPRINT/ --- files_x\StartEnd\file_start2.txt ---
$$-> CSYS / 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,  $
            0.0000000000, 1.0000000000, 0.0000000000, 0.0000000000,  $
            0.0000000000, 0.0000000000, 1.0000000000, 0.0000000000
MULTAX/ ON
MODE / MILL
SPINDL / MAIN
PPRINT/ --- files_x\ToolChange\file_TC_start.txt ---
PPRINT/ --- files_x\Comment\toolChange_comment.txt --- 
PPRINT/ - T1 R    4.00000 L   16.00000     0.00000 WAAM TORCH2  
PPRINT/ OPERATION 3 
PPRINT/  
PPRINT/ Mold Repair 
PPRINT/ T1 Additive Manufacturing 
PPRINT/ --- 
PPRINT/ --- files_x\Additive\additive_toolchange.txt ---
$$ LASER TOOL CHANGE
PPRINT/ ---
PPRINT/ --- files_x\ToolChange\loadTool.txt ---
LOADTL/1
PPRINT/ --- files_x\ToolChange\file_TC_end.txt ---
PPRINT/ --- files_x\Comment\job_start.txt ---
PPRINT/ OPERATION 3 
PPRINT/  
PPRINT/ Mold Repair 
PPRINT/ T1 Additive Manufacturing 
CALSUB/START_JOB
SEQUENCE/ BEGIN,toolpath
PPRINT/ --- 
FEDRAT/ MMPM, 450
PPRINT/ --- files_x\Additive\additive_technology_peripheral.txt ---
$$ ADDITIVE TECHNOLOGY PHERIPHERAL
SPINDL/OFF
PPRINT/ ---
PPRINT/ --- files_x\Additive\additive_layer.txt ---
$$ Layer: 1
PPRINT/ --- 
PPRINT/ --- files_x\Jobfiles\5X_begin.txt ---
RAPID
GOTO /  -51.86086,  -37.57609,  130.00000, 0.000000, 0.000000, 1.000000 ,    0.00000 ,    0.00000
RAPID
GOTO /  -51.86086,  -37.57609,  112.35300, 0.000000, 0.000000, 1.000000 ,    0.00000 ,    0.00000
PPRINT/ --- files_x\Additive\additiveDevice_on.txt ---
CALSUB/START_DEPO
PPRINT/ ---
FEDRAT/ MMPM, 450
GOTO /  -45.36639,  -45.18014,  112.35300, 0.000000, 0.000000, 1.000000 ,    0.00000 ,    0.00000
PPRINT/ --- files_x\Additive\additiveDevice_off.txt ---
CALSUB/STOP_DEPO
PPRINT/ ---
PPRINT/ --- files_x\Additive\additive_technology_peripheral.txt ---
$$ ADDITIVE TECHNOLOGY PHERIPHERAL
SPINDL/OFF
PPRINT/ ---
RAPID
GOTO /  -45.36639,  -45.18014,  130.00000, 0.000000, 0.000000, 1.000000 ,    0.00000 ,    0.00000
PPRINT/ --- files_x\Additive\additive_layer.txt ---
$$ Layer: 1
PPRINT/ --- 
RAPID
GOTO /  -29.17694,   45.00831,  130.00000, 0.000000, 0.000000, 1.000000 ,    0.00000 ,    0.00000
RAPID
GOTO /  -29.17694,   45.00831,  112.35300, 0.000000, 0.000000, 1.000000 ,    0.00000 ,    0.00000
PPRINT/ --- files_x\Additive\additiveDevice_on.txt ---
CALSUB/START_DEPO
PPRINT/ ---
FEDRAT/ MMPM, 450
GOTO /  -22.68247,   37.40425,  112.35300, 0.000000, 0.000000, 1.000000 ,    0.00000 ,    0.00000
PPRINT/ --- files_x\Additive\additiveDevice_off.txt ---
CALSUB/STOP_DEPO
PPRINT/ ---
RAPID
GOTO /  -22.68247,   37.40425,  130.00000, 0.000000, 0.000000, 1.000000 ,    0.00000 ,    0.00000
PPRINT/ --- files_x\Jobfiles\5X_end.txt ---
PPRINT/ --- files_x\Comment\job_end.txt ---
SEQUENCE/ END
CALSUB/END_JOB
PPRINT/ --- files_x\StartEnd\file_end.txt ---
SPINDL/OFF
CALSUB/END_PROG
FINI
