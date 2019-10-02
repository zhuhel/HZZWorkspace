import ROOT 
from numpy import linspace 

# Access ROOT file 
f = ROOT.TFile( '/home/goblirsc/Code/H4l/Workspaces_CI/build/combined.root' ) 

# Access the RooWorkspace 
w = f.Get( 'combined' )

# Print the workspace contents 
w.Print() 

# Retrieve the simulatenous PDF 
pdf = w.pdf( 'simPdf' ) 

# Retrieve the Asimov dataset 
asimov = w.data( 'MC' ) 

# Retrieve m4l
m4l = w.obj( 'gaussThing' ) 

# Retrieve mH 
mainObs = w.obj( 'mean' ) 
mainObs.Print() 

##::: Draw the PDF 
#----------------------
# Create a canvas 
canvas = ROOT.TCanvas( "combined_template" ) 

# # Create a RooFit frame from the m4l observable 
# m4lframe = m4l.frame() 

# color = 39 
# for value in linspace( -5, 5, 20 ): 
#   # Set the value of mH and make it constant  
#   mainObs.setVal( value )
#   mainObs.setConstant( True ) 
  
#   # Plot the PDF 
#   pdf.plotOn( m4lframe, ROOT.RooFit.ProjWData( asimov ), ROOT.RooFit.LineColor(   color + 1 ) )
#   color = color + 1

# # Draw the RooFit frame 
# m4lframe.Draw() 

# # Save the canvas 
# canvas.SaveAs( ".pdf" )

# Fit the asimov dataset with the PDF 
mainObs.setConstant( False ) 
mainObs.setVal( 1 )
pdf.fitTo( asimov, ROOT.RooFit.SumW2Error( ROOT.kTRUE ) )
mainObs.Print()
