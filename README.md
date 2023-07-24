Project README
This project consists of a series of scripts for processing and extracting specific information from different CSV files and writing the selected information into new text files. Here is a brief description of the scripts:

FE Mesh Nodes Extraction: This script reads the 'FE Mesh Nodes.csv' file, extracts the 'No.', 'X [mm]', 'Y [mm]', and 'Z [mm]' columns, and writes this information into a new text file called 'Nodes.txt'.

Cross-Sections Properties Extraction: This script reads the '1.13 Cross-Sections.csv' file, extracts the 'Torsion J', 'Bending Iy', 'Bending Iz', and 'Axial A' columns, and writes this information into a new text file called 'Members.txt'.

Cross-Sections Names Extraction: This script reads the '1.13 Cross-Sections.csv' file again, this time extracts the 'Description [mm]' column, and writes the names into a new text file called 'Names.txt'.

FE Mesh Cells Extraction: This script reads the 'FE Mesh Cells .csv' file, extracts the 'No.', 'No..1', '1', '2', '3', and '4' columns, and writes this information into a new text file called 'Mesh_Cells.txt'.

Surfaces Extraction: This script reads the '1.4 Surfaces.csv' file, extracts the 'No.', 'No..1', and 'd [mm]' columns, and writes this information into a new text file called 'Surface.txt'.

Element ID Extraction: This script reads the '1.17 Members.csv' file, extracts the 'No.' and 'Start' columns, and writes this information into a new text file called 'Elementid.txt'.

Lines Extraction: This script reads the '1.2 Lines.csv' file, splits the 'Nodes No.' column into 'Start Node' and 'End Node' columns, and writes this information along with the 'No.' column into a new text file called 'Lines.txt'.

These scripts are run in Python, utilizing the pandas library for data manipulation. The generated text files are used in a 3D model for further computations.

File Headers Information
In this project, specific data are extracted from CSV files and saved into text files. These text files are used in a 3D model for further computations. Here are the headers for each text file created:

Nodes.txt

Column 1: Node Number (No.)
Column 2: X Coordinate (X [mm])
Column 3: Y Coordinate (Y [mm])
Column 4: Z Coordinate (Z [mm])
Members.txt

Column 1: Torsion Constant (Torsion J)
Column 2: Bending Moment of Inertia along Y-axis (Bending Iy)
Column 3: Bending Moment of Inertia along Z-axis (Bending Iz)
Column 4: Cross-sectional Area (Axial A)
Names.txt

Column 1: Cross-Section Description (Description [mm])
Mesh_Cells.txt

Column 1: Mesh Cell Number 1 (No.)
Column 2: Mesh Cell Number 2 (No..1)
Column 3: Node 1 (1)
Column 4: Node 2 (2)
Column 5: Node 3 (3)
Column 6: Node 4 (4)
Surface.txt

Column 1: Surface Number 1 (No.)
Column 2: Surface Number 2 (No..1)
Column 3: Thickness (d [mm])
Elementid.txt

Column 1: Member Number (No.)
Column 2: Start Node Number (Start)
Lines.txt

Column 1: Line Number (No.)
Column 2: Start Node (Start Node)
Column 3: End Node (End Node)

