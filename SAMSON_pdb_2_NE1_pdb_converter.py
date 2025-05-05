#Here
input_file = open('flyingH.pdb', 'r') 
igot = input_file.readlines()
f = open('flyingHNE1.pdb', "a") 
#this is the pdb version format
#COLUMNS      DATA TYPE        FIELD      DEFINITION
#------------------------------------------------------
 #1 -  6      Record name      "ATOM    "
 #7 - 11      Integer          serial     Atom serial number.
#13 - 16      Atom             name       Atom name.
#17           Character        altLoc     Alternate location indicator.
#18 - 20      Residue name     resName    Residue name.
#22           Character        chainID    Chain identifier.
#23 - 26      Integer          resSeq     Residue sequence number.
#27           AChar            iCode      Code for insertion of residues.
#31 - 38      Real(8.3)        x          Orthogonal coordinates for X in
                                         #Angstroms
#39 - 46      Real(8.3)        y          Orthogonal coordinates for Y in
                                         #Angstroms
#47 - 54      Real(8.3)        z          Orthogonal coordinates for Z in
                                         #Angstroms
#55 - 60      Real(6.2)        occupancy  Occupancy.
#61 - 66      Real(6.2)        tempFactor Temperature factor.
#77 - 78      LString(2)       element    Element symbol, right-justified.
#79 - 80      LString(2)       charge     Charge on the atom.
#HETATM    1  C       A   1       7.629  -1.386 -23.948  1.00  0.00           C
for line_sam in igot:
    if line_sam.find("HETATM ") > -1:        
        SAM_atoms =  line_sam.split()
        tag = "HETATM" #SAM_atoms[0]
        index = SAM_atoms[1]
        element = SAM_atoms[8]
        xcor = SAM_atoms[5]
        ycor = SAM_atoms[6]
        zcor = SAM_atoms[7]
        
        record_name_len = 6
        record_name_pad = record_name_len - len(tag)
        
        serial_len = 5
        index_pad = serial_len - len(index) 
        
        atom_name_len = 3
        atom_name_pad = atom_name_len - len(element)
        
        character = 1
        residue_name = 3
        
        chainID_len =   4
        chainID = "A" 
        chainID_pad = chainID_len - len(chainID)
        
        residue_sequence_num_len = 4
        residue_sequence_num = "1"
        residue_seqeunce_num_pad = residue_sequence_num_len- len(residue_sequence_num)
        
        iCode = 4
        
        x_len = 8
        x_len_pad = x_len - len(xcor)
        
        y_len = 8
        y_len_pad = y_len - len(ycor)        
        
        z_len = 8
        z_len_pad = z_len - len(zcor)        
       
        
        ##all text end at the rightside of its column allowance, so there wouldn't be padding on the right side, the left side padding just descreases and the string gets longer.
        pdb_line = record_name_pad *" " + tag +index_pad*" " +index + atom_name_pad*" " + element + character*" " + residue_name* " " + chainID_pad*" " + chainID + residue_seqeunce_num_pad*" " + residue_sequence_num + iCode*" "+ x_len_pad*" "+ xcor + y_len_pad*" "+ ycor + z_len_pad*" "+ zcor+"  1.00  0.00           "+ element       +"\n"
        print(pdb_line, len(pdb_line))
        f.write(pdb_line)
        
    #this needs a more elegant solution for larger atom counts
    if line_sam.find("HETATM1") > -1:
        print(line_sam)
        SAM_atoms =  line_sam.split()
        tag = SAM_atoms[0]
        #index = SAM_atoms[1]
        element = SAM_atoms[7]
        xcor = SAM_atoms[4]
        ycor = SAM_atoms[5]
        zcor = SAM_atoms[6]
        
        record_name_len = 6
        record_name_pad = record_name_len - len(tag)
        
        serial_len = 5
        index_pad = serial_len - len(index) 
        
        atom_name_len = 3
        atom_name_pad = atom_name_len - len(element)
        
        character = 1
        residue_name = 3
        
        chainID_len =   4
        chainID = "A" 
        chainID_pad = chainID_len - len(chainID)
        
        residue_sequence_num_len = 4
        residue_sequence_num = "1"
        residue_seqeunce_num_pad = residue_sequence_num_len- len(residue_sequence_num)
        
        iCode = 4
        
        x_len = 8
        x_len_pad = x_len - len(xcor)
        
        y_len = 8
        y_len_pad = y_len - len(ycor)        
        
        z_len = 8
        z_len_pad = z_len - len(zcor)        
       
        
        ##all text end at the rightside of its column allowance, so there wouldn't be padding on the right side, the left side padding just descreases and the string gets longer.
        pdb_line = record_name_pad *" " + tag + atom_name_pad*" " + element + character*" " + residue_name* " " + chainID_pad*" " + chainID + residue_seqeunce_num_pad*" " + residue_sequence_num + iCode*" "+ x_len_pad*" "+ xcor + y_len_pad*" "+ ycor + z_len_pad*" "+ zcor+"  1.00  0.00           "+ element       +"\n"
        print(pdb_line, len(pdb_line))
        f.write(pdb_line)        
    
        
    
    
    
    if line_sam.find("CONECT") > -1: 
        f.write(line_sam)












input_file.close()
f.close()