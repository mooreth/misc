def get_range_of_velocities(path):
    """ Returns the range of velocity norm """
    minV = maxV = SBVelocity3().norm()
    if path:
        if path.hasVelocityData():
            velocityData = path.getVelocityData()            
            if len(velocityData):
                if len(velocityData[0]):
                    minV = maxV = velocityData[0][0].norm()
                # go through frames
                for velocities in velocityData:
                    # go through velocities per frame
                    for velocity in velocities:
                        vn = velocity.norm() #.norm returns the normal of the velocity vectors
                        if vn > maxV: maxV = vn
     
                        elif vn < minV: minV = vn
 
    
    return [minV, maxV] #for the c60 collision at 2 these equal 0 pm.fs^-1 & 105.699 pm.fs^-1
    
    

# get path
scaling_factor = 10 #If all your numbers are real small e.g. 0.01, 0.02, this can make the final captures more colorful
pathIndexer = SAMSON.getNodes('node.type path')

if len(pathIndexer):
    # get the path from the user
    ret, capture_filepath = SAMSON.getPathFromUser(dialogTitle = 'Choose folder to save captures...')
else:
    ret = False

if ret and len(pathIndexer):
    # get the 1st path
    path = pathIndexer[0]
    if path.hasVelocityData():
        # get all atoms in the path
        atomIndexer = path.getAtomIndexer()

        minVelNorm, maxVelNorm = get_range_of_velocities(path)
        minVelNorm = SBVelocity3().norm() # set the min value to 0 ####I assume is the min velocity is not 0, this sets to the bottom of the velocity scale i.e. 0s it

        # choose a default color palette
        palette = SBPaletteDefaultPalette.discreteTab20c

        # make the operation undoable
        with SAMSON.holding("Colorize atoms"):
            for step in range(0, len(path)):

                path.currentStep = step

                # colorize each atom according to its coordinate
                for atom in atomIndexer:
                    velocity = path.getVelocity(step, atom)
                    #print(velocity[0].value)
                    #print("atom velocity", velocity) #this is the velocity in the trajectory file, listed by atom index number
                    velNorm = velocity.norm()
                    #print("normalized velocity", velNorm)
                    # get the color 'intensity' value, which should be in the [0, 1] range
                    
                    h = ((velNorm - minVelNorm) / (maxVelNorm - minVelNorm)).value * scaling_factor
                    
                    # get RGB color from the palette
                    color = palette.getColor(h)
                    # set the color of the node
                    atom.setColor(color)
                
                SAMSON.captureViewportToFile(
                    filename = f"{capture_filepath}/capture-{step}.png",
                    width = 1200, height = 800,
                    transparentBackground = False, usePathTracing = False, showProgressBar = False)