#include "TSystem.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "HZZWorkspace/PlotHelp.h"
#include "Rtypes.h"
#include <map>
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TLine.h"
#include "TLatex.h"
#include "TMath.h"
#include "TROOT.h"

#include <iostream>
#include <fstream>
#include <stdio.h>

#include "TXMLEngine.h"

//------------------------------------------------------------------------------
// Building acceptance polynomial for given minitree input and configuration (xml format)
//
// Hardcoded configurations commented out at bottom of the file (temporary!)
//------------------------------------------------------------------------------
int main( int argc, char *argv[] ){

  // Get the name of the configuration file (in xml format) from the arguments
  if( argc < 1 || std::string( argv[ 1 ] ).find( ".xml" ) != std::string::npos ){
    std::cout << "Configuration file missing. Please see examples/Tools/ac_massmeasurement_config.xml for an example" << std::endl;
    std::cout << "More information available here: https://gitlab.cern.ch/HZZ/HZZSoftware/HZZWorkspace/wikis/Executables/makeACpolynomials" << std::endl;
    std::cerr << "Please provide a configuration in xml format" << std::endl;
  }
  std::string config = argv[ 1 ];

  if( argc > 2 ){ std::cout << "I don't know what to do with more than 1 passed argument" << std::endl; }

  // Declare string for minitree directory
  std::string minitreeDir;

  // Declare a map to store the mH upper and lower boundaries
  std::map< std::string, float > mH{
    { "min", 0 },
    { "max", 0 }
  };

  // Declare a map for the plotting boundaries
  std::map< std::string, float > plotBoundaries{
    { "min", 0 },
    { "max", 0 },
    { "minY", 0.05 }, // Set to -1 to use auto
    { "maxY", 0.2 } // Set to -1 to use auto
  };

  // Establish the output plot images extension
  std::string imgExt = "eps";

  // Declare a vector of strings to store the width information
  std::vector< std::string > widths;

  // Declare vectors of strings for the production mode information
  std::vector< std::string >* prod = new std::vector< std::string >();
  std::vector< std::string >* prodlabel = new std::vector< std::string >();

  // Create a map to associate a string with each vector
  std::map< std::string, std::vector< std::string >* > ProductionModes{
    { "name", prod },
    { "label", prodlabel }
  };

  // Declare vectors of strings for the category information
  std::vector< std::string >* categories = new std::vector< std::string >();
  std::vector< std::string >* catlabels = new std::vector< std::string >();
  std::vector< std::string >* cuts = new std::vector< std::string >();
  std::vector< std::string >* colors = new std::vector< std::string >();
  std::vector< std::string >* markers = new std::vector< std::string >();

  // Create a map to associate a string with each vector
  std::map< std::string, std::vector< std::string >* > Categories{
    { "name", categories },
    { "label", catlabels },
    { "selection", cuts },
    { "color", colors },
    { "marker", markers }
  };

  // Open an XML Engine to read the configuration file
  TXMLEngine xml;

  // Open the XML document
  XMLDocPointer_t xmldoc = xml.ParseFile( config.c_str() );

  // Access the main node of the document
  XMLNodePointer_t mainnode = xml.DocGetRootElement( xmldoc );

  // Retrieve the first child node
  XMLNodePointer_t childnode = xml.GetChild( mainnode );

  // Loop on the child nodes to populate the configuration
  while( childnode != 0 ){

    // Get the name of the node
    std::string nodename = xml.GetNodeName( childnode );

    // Populate the minitree directory
    if( nodename.find( "minitreedir" ) != std::string::npos ){
      minitreeDir = std::string( xml.GetNodeContent( childnode ) );
      // Add the last directory / if it's missing
      if( minitreeDir.back() != '/' ) minitreeDir += "/";
      // Report minitree directory to  commandline  
      std::cout << "Looking for files in " << minitreeDir << std::endl;
    }

    // Populate widths vector
    if( nodename.find( "width" ) !=  std::string::npos ){
      widths.push_back( std::string( xml.GetNodeContent( childnode ) ) );
    }

    // Access grandchild nodes -- "details"
    XMLNodePointer_t details = xml.GetChild( childnode );

    // If there are details, process them
    while( details != 0 ){
      // Identify this subnode
      std::string property = xml.GetNodeName( details );
      std::string value = std::string( xml.GetNodeContent( details ) );

      // Store mH properties
      if( nodename == "mH" && mH.find( property ) != mH.end() ){
        mH[ property ] = std::atof( value.c_str() );
        std::cout << Form( "Found %s:%s with value %s", nodename.c_str(), property.c_str(), value.c_str() ) << std::endl;
      }

      // Store production mode properties
      if( nodename == "productionmode" && ProductionModes.find( property ) != ProductionModes.end() ){
        ProductionModes[ property ]->push_back( value );
        std::cout << Form( "Found %s:%s with value %s", nodename.c_str(), property.c_str(), value.c_str() ) << std::endl;
      }

      // Store category properties
      if( nodename == "category" && Categories.find( property ) != Categories.end() ){
        Categories[ property ]->push_back( value );
        std::cout << Form( "Found %s:%s with value %s", nodename.c_str(), property.c_str(), value.c_str() ) << std::endl;
      }

      // Advance to the next details grandchild
      details = xml.GetNext( details );
    }
    // Get the next child node
    childnode = xml.GetNext( childnode );
  }

  // Release the XML from memory
  xml.FreeDoc( xmldoc );

  // Convert strings to floats where necessary
  float minMH = mH[ "min" ];
  float maxMH = mH[ "max" ];
  float minPlot = plotBoundaries[ "min" ];
  float maxPlot = plotBoundaries[ "max" ];
  float minPlotY = plotBoundaries[ "minY" ];
  float maxPlotY = plotBoundaries[ "maxY" ];

  std::cout << minMH << std::endl;
  std::vector< int > color;
  for( auto i : *colors ){
    color.push_back( atoi( i.c_str() ) );
  }
  std::vector< int > marker;
  for( auto i : *markers ){
    marker.push_back( atoi( i.c_str() ) );
  }

  // float lumi=1;
  // bool storedLumi=false;

  //::: SET UP RESULTS FILES
  // Establish file names
  std::string outTFileName = "ACplots.root";
  std::string outFileName = "polyNorm.txt";
  std::string outFileSIName = "polyNormSI.txt";
  std::string outFileSysName = "polySys.txt";

  // Create a directory for the plots directory
  gSystem->Exec( "mkdir -p $PWD/plots/" );

  // Create and open output ROOT file
  TFile* outTfile = new TFile( Form( "plots/%s", outTFileName.c_str() ), "RECREATE" );

  // Create output text files
  std::ofstream outfile;
  std::ofstream outfileSI;
  std::ofstream outfileSys;

  // Open output text files
  outfile.open( outFileName.c_str(), std::ios::out );
  outfileSI.open( outFileSIName.c_str(), std::ios::out );
  outfileSys.open( outFileSysName.c_str(), std::ios::out );

  //
  // LOOP OVER PRODUCTION MODES AND WIDTHS77
  //
  // Create a canvas for the acceptance plot
  TCanvas can( "can", "", 600, 600 );
  can.SetGrid( 2,2 );
  can.SetLeftMargin( 0.15 );

  // Loop on production mode
  for( unsigned int p( 0 ); p < prod->size(); ++p ){
    // Loop on width
    for( unsigned int w( 0 ); w < widths.size(); ++w ){

      // Production mode and width assigned
      std::cout << "\nprod: " << prod->at( p ) << " width: " << widths[ w ] << std::endl;

      // Based on the width keyword, decide whether taus are included
      bool hasTau = true;
      if( widths[ w ] == "NW" ) hasTau = true;
      if( widths[ w ] == "ZZ4lep" ) hasTau = false;

      //
      // BUILD MAP OF (MASS POINT VS FILE)
      //
      std::map< float, std::string > fileMap;
      std::string str;
      const char* entry;

      // Open the minitree directory
      void* dirp = gSystem->OpenDirectory( minitreeDir.c_str() );

      // Loop on all files in the minitree directory
      while( ( entry = ( char* )gSystem->GetDirEntry( dirp ) ) ){
        // Name of minitree
        str = entry;
        size_t pos = str.find( prod->at( p ).c_str() );
        size_t pos2 = str.find( widths[ w ].c_str() );

        if( pos != std::string::npos && pos2 != std::string::npos ){
          pos += prod->at( p ).size()-1;

          // Identify the masspoint
          std::string massString = "";
          while( !isdigit( str[ ++pos ] ) ){}
          while( isdigit( str[ pos ] ) ){ massString += str[ pos++ ]; }
          if( str.substr( pos, 2 ) == "p5" ) massString += ".5";
          if( str.substr( pos, 3 ) == "p25" ) massString += ".25";
          float mass = atof( massString.c_str() );

          // Mass point lies outside the range of interest
          if( mass < minMH || mass > maxMH ) continue;

          // Report file and masspoint to the command line
          std::cout << "From file: " << str << " --> extracted mass = " << mass << std::endl;

          // Check against the tau condition
          if( hasTau && str.find( "noTau" ) != std::string::npos ) continue;
          if( !hasTau && str.find( "noTau" ) == std::string::npos ) continue;

          // Store the minitree in the file map
          fileMap[ mass ] = minitreeDir + str;
        }
      }
      // Report map built and give the number of entries
      std::cout << "Built a map of file names vs mass. Size = " << fileMap.size() << std::endl;

      // No files found? Exit loop
      if( fileMap.size() == 0 ) continue;

      //BUILD A MAP OF (MASS POINT VS ACC)
      std::map< std::string, std::vector< float > > norm, norm_error;
      std::vector< float > zero, masses;

      // Loop on the selected minitrees
      for (std::map< float, std::string >::iterator it = fileMap.begin(); it != fileMap.end(); ++it ){

        std::cout << "For mass = " << it->first << ", using file " << it->second << std::endl;

        // Open a new ROOT file
        TFile* file = new TFile( it->second.c_str(), "READ" );

        // Access the minitree
        TTree* tree = ( TTree* )file->Get( "tree_incl_all" );

        masses.push_back( ( *it ).first);
        zero.push_back( 0. );

        // Loop on the categories
        for( unsigned int c( 0 ); c < categories->size(); ++c ){
          // Createa a string for the current category
          std::string cat = categories->at( c );
          // Create a histogram with a single bin to serve as a counter
          TH1F* h = new TH1F( "h", "h", 1, -9999, 9999 );
          // Populate it using TTree::Draw(). Select events matching this category
          tree->Draw( "m4l_constrained_HM>>h", Form( "%s*(9./4.)*(weight/(w_lumi*w_xs*w_br))", cuts->at( c ).c_str() ) );
          // Store the count and uncertainty
          norm[ cat ].push_back( h->GetBinContent( 1 ) );
          norm_error[ cat ].push_back( h->GetBinError( 1 ) );
          // Remove the histogram
          delete h;
        }
      }

      // Report that the acceptance map is complete
      std::cout << "Built a map of acceptance vs mass. Size = " << masses.size() << std::endl;

      //::: Report the acceptances to the commandline
      std::cout << "MC acceptances for " << prod->at( p ) << ":" << std::endl;
      // Loop on  the categories
      for( unsigned int c( 0 ); c < categories->size(); ++c ){
        // Name the category
        std::cout << "\t" << categories->at( c ) << std::endl;
        // Loop on the mass points
        for( unsigned int m( 0 ); m < masses.size(); ++m ){
          // Print the acceptance and its uncertainty
          std::cout << "\t\t" << masses[ m ] << ":\t" << norm[ categories->at( c ) ][ m ] << "\t+-\t" << norm_error[ categories->at( c ) ][ m ] << std::endl;
        }
      }

      // Create a character string to identify the polynomial order
      const char* pol="pol0";
      if( masses.size() > 1 ) pol = "pol1";
      if( masses.size() > 2 ) pol = "pol2";
      //if (masses.size()>3) pol="pol3";//FIXME
      //if (masses.size()>4) pol="pol4";
      //if (masses.size()>5) pol="pol5";
      //if (masses.size()>5) pol="pol6";

      //:::: BUILD A GRAPH OF MASS POINT VS NORMS
      // Create a new multigraph
      TMultiGraph* mg = new TMultiGraph();
      // Access the canvas' draw space
      can.cd();
      // Loop on the categories
      for( unsigned int c( 0 ); c < categories->size(); ++c ){
        // Identify the current category
        std::string cat = categories->at( c );

        // Create a new graph with uncertainties where x is the mass point and y is the acceptance wit uncertainty
        TGraphErrors* graph = new TGraphErrors( masses.size(), &( masses[ 0 ] ), &( norm[ cat ][ 0 ] ), &( zero[ 0 ] ), &( norm_error[ cat ][ 0 ] ) );
        // Fit the graph with the polynomial. Perform a quiet fit with minimum printing and return the result
        TFitResultPtr fit = graph->Fit( pol, "qs" );
        // Add this category's graph to the multi-graph
        mg->Add( graph );

        // Enter the output ROOT file
        outTfile->cd();
        // Give a name to the graph so that it is uniquely written
        std::string graphName = Form( "AC_%s_%s_%s", prod->at( p ).c_str(), widths[ w ].c_str(), cat.c_str() );
        // Save the graph in the ROOT file
        graph->Clone( graphName.c_str() )->Write();

        //::: ESTIMATE SYSTEMATICS
        // Extract the fitted TF1 from the graph
        TF1* polfit = graph->GetFunction( pol );

        // Create a vector to hold the deviations
        std::vector< float > dev;
        // Loop on the mass points
        for( unsigned int m( 0 ); m < masses.size(); ++m ){
          // Calculate the difference between the fitted function at the masspoint and its value from Monte Carlo
          dev.push_back( polfit->Eval( masses[ m ] ) - norm[ cat ][ m ] );
        }
        // Calculate the RMS, standardized with respect to the number of mass points, of the deviations
        float rms = TMath::RMS( dev.size(), &( dev[ 0 ] ) )/sqrt( masses.size()-1 );

        // Clear out the vector
        dev.clear();

        //  Save the deviation to the systematics file
        for( unsigned int m( 0 ); m < masses.size(); ++m ){
          float nom = polfit->Eval( masses[ m ] );
          outfileSys << prod->at( p ) << "  " << cat << "  " << masses[ m ] << "  " << Form( "%.8f", ( nom-rms )/ nom ) << "   " << Form( "%.8f", ( nom+rms )/nom ) << std::endl;
        }

        // Printout for conf note table
        std::cout << "prod " << prod->at( p ) << " in category " << cat << " acceptance @125GeV = " << 100*polfit->Eval( 125 ) << "%,    @126GeV = " << 100*polfit->Eval( 126 ) << "%" << std::endl;

        // Save the graph as an image
        // Enter the canvas
        can.cd();
        // graph's general properties
        graph->SetMinimum( 0 );
        graph->SetTitle( "" );
        // Function properties
        graph->GetFunction( pol )->SetLineColor( color[ c ] );
        // x axis
        graph->GetXaxis()->SetTitle( "m_{H} [GeV]" );
        graph->GetXaxis()->SetRangeUser( masses[ 0 ], masses[ masses.size()-1 ] );
        graph->GetXaxis()->CenterTitle();
        graph->GetXaxis()->SetTitleFont( 42 );
        graph->GetXaxis()->SetTitleSize( 0.04 );
        graph->GetXaxis()->SetLabelSize( 0.035 );
        graph->GetXaxis()->SetRangeUser( minPlot, maxPlot );
        //  y axis
        graph->GetYaxis()->SetTitle( Form( "%s  acceptance", prodlabel->at( p ).c_str() ) );
        graph->GetYaxis()->SetTitleFont( 42 );
        graph->GetYaxis()->SetTitleOffset( 1.4 );
        graph->GetYaxis()->SetTitleSize( 0.04 );
        graph->GetYaxis()->SetLabelSize( 0.035 );
        if( minPlotY != -1 ) graph->SetMinimum( minPlotY );
        if( maxPlotY != -1 ) graph->SetMaximum( maxPlotY );
        // marker style
        graph->SetMarkerStyle( marker[ c ] );
        graph->SetMarkerSize( 0.75 );
        graph->SetMarkerColor( color[ c ] );
        // line style
        graph->SetLineColor( color[ c ] );
        // Draw the graph with axes and points
        graph->Draw( "ap" );

        // Create a legend manually
        //    Create a line matched to the category's TGraph
        TLine* line = new TLine( 0.65, 0.75, 0.68, 0.75 );
        line->SetNDC( true );
        line->SetLineColor( color[ c ] );
        line->SetLineWidth( 3 );
        line->Draw( "same" );
        //    Label the line with the category label
        myText( 0.57, 0.74, kBlack, catlabels->at( c ).c_str(), 0.04 );
        // Add the ATLAS designation
        //ATLASLabel( 0.18, 0.84, "Internal", kBlack );

        // Create a directory for the plots output
        gSystem->Exec( Form( "mkdir -p $PWD/plots/acceptance/%s/%s/", prod->at( p ).c_str(), widths[ w ].c_str() ) );
        // Save the canvas as an image
        can.Print( Form( "$PWD/plots/acceptance/%s/%s/acceptance_%s_%s_%s_%s.%s", prod->at( p ).c_str(), widths[ w ].c_str(), prod->at( p ).c_str(), widths[ w ].c_str(), cat.c_str(), pol, imgExt.c_str() ) );

        //::: STORE POLY PARAMS
        std::vector< double > polyParameters = fit->Parameters();
        // Create header lines within the output text files
        outfile << Form( "[%s %s %s]\n", prod->at( p ).c_str(), widths[ w ].c_str(), cat.c_str() );
        outfileSI << Form( "[%s %s %s]\n", prod->at( p ).c_str(), widths[ w ].c_str(), cat.c_str() );
        // Create the separators
        std::string outline = "";
        std::string outlineSI = "& ";
        // Add the parameters to the line
        for( unsigned int r( 0 ); r < polyParameters.size(); ++r ){
          outline = outline + " " + Form( "%.20f", polyParameters[ r ] );
          outlineSI = outlineSI + Form( "%.3e", polyParameters[ r ] ) + " & ";
        }
        // Add new line characters
        outlineSI += "\n";
        outline += "\n";
        // Record the lines
        outfile << outline;
        outfileSI << outlineSI;
      }
      // Draw the multigraph  plot
      mg->Draw("ap");
      // x axis
      mg->GetXaxis()->SetRangeUser( masses[ 0 ], masses[ masses.size()-1 ] );
      mg->GetXaxis()->SetTitle( "m_{H} [GeV]" ); //m_S for high mass //FIXME
      mg->GetXaxis()->CenterTitle();
      mg->GetXaxis()->SetTitleFont( 42 );
      mg->GetXaxis()->SetTitleSize( 0.04 );
      mg->GetXaxis()->SetLabelSize( 0.035 );
      mg->GetXaxis()->SetRangeUser( minPlot, maxPlot );
      // y axis
      if( minPlotY != -1 ) mg->SetMinimum( minPlotY );
      if( maxPlotY != -1 ) mg->SetMaximum( maxPlotY );
      mg->GetYaxis()->SetTitle(Form( "%s  acceptance", prodlabel->at( p ).c_str() ) );
      mg->GetYaxis()->SetTitleFont( 42 );
      mg->GetYaxis()->SetTitleOffset( 1.4 );
      mg->GetYaxis()->SetTitleSize( 0.04 );
      mg->GetYaxis()->SetLabelSize( 0.035 );

      // Add labels
      float x( 0.25 ), y( 0.85 );
      //   ATLAS
      ATLASLabel( x, y, "Internal", kBlack );
      //   Categories with their lines
      for (unsigned int c( 0 ); c < categories->size(); ++c )
          myLineText( x, y-0.05-0.04*c, color[ c ], kSolid, 2, catlabels->at( c ).c_str(), 0.035 );

      // Update the canvas so that changes take effect
      can.Update();
      // Save the canvas as an image
      can.Print( Form( "$PWD/plots/acceptance/%s/%s/acceptance_%s_%s_%s.%s", prod->at( p ).c_str(), widths[ w ].c_str(), prod->at( p ).c_str(), widths[ w ].c_str(), pol, imgExt.c_str() ) );

    } //end loop over widths

  } //end loop over prods

  // Close all output files
  outTfile->Close();
  outfile.close();
  outfileSI.close();
  outfileSys.close();
}
/*
  float minMH = 120;  //Range to user MC mass points into fit //FIXME
  float maxMH = 130;
  float minPlot = 122; //Range to plot fits //FIXME MELA uses 300, LWA uses 400,  NWA uses 200
  float maxPlot = 128;
  float minPlotY=0.05; //set to -1 to use auto
  float maxPlotY=0.2; //set to -1 to use auto

  // float lumi=1;
  // bool storedLumi=false;

  std::vector<std::string> widths;
  //widths.push_back("NW");
  widths.push_back("ZZ4lep");
  //widths.push_back("w5");
  //widths.push_back("w10");
  //widths.push_back("w15");

  std::vector<std::string> prod;
  prod->push_back("ggH");
  //prod->push_back("HZZ4l");
  //prod->push_back("VBFH");

  std::vector<std::string> prodlabel;
  prodlabel.push_back("ggF");
  //prodlabel.push_back("ggF (15\%)");
  //prodlabel.push_back("VBF");

  std::vector<std::string> categories;
  std::vector<std::string> cuts;
  std::vector<std::string> catlabels;
  std::vector<int> color;
  std::vector<int> marker;


  //Mass Msmt
  categories->push_back("ggF_4mu"); color.push_back(kViolet); marker.push_back(24); catlabels.push_back("4#mu");
  categories->push_back("ggF_4e"); color.push_back(kOrange); marker.push_back(25); catlabels.push_back("4e");
  categories->push_back("ggF_2mu2e"); color.push_back(kGreen+2); marker.push_back(26); catlabels.push_back("2#mu2e");
  categories->push_back("ggF_2e2mu"); color.push_back(kCyan); marker.push_back(27); catlabels.push_back("2e2#mu");
  cuts.push_back("(pass_vtx4lCut==1 && 110<m4l_constrained&&m4l_constrained<135&&event_type==0)");
  cuts.push_back("");
  cuts.push_back("(pass_vtx4lCut==1 && 110<m4l_constrained&&m4l_constrained<135&&event_type==2)");
  cuts.push_back("(pass_vtx4lCut==1 && 110<m4l_constrained&&m4l_constrained<135&&event_type==3)");
  //cuts.push_back("(pass_vtx4lCut==1 && 110<m4l_fsr&&m4l_fsr<135&&event_type==0)");
  //cuts.push_back("(pass_vtx4lCut==1 && 110<m4l_fsr&&m4l_fsr<135&&event_type==1)");
  //cuts.push_back("(pass_vtx4lCut==1 && 110<m4l_fsr&&m4l_fsr<135&&event_type==2)");
  //cuts.push_back("(pass_vtx4lCut==1 && 110<m4l_fsr&&m4l_fsr<135&&event_type==3)");

  //NWA
  categories->push_back("ggF_4mu"); color.push_back(kViolet); marker.push_back(24); catlabels.push_back("4#mu ggF-like");
  categories->push_back("ggF_4e"); color.push_back(kOrange); marker.push_back(25); catlabels.push_back("4e ggF-like");
  categories->push_back("ggF_2mu2e"); color.push_back(kGreen+2); marker.push_back(26); catlabels.push_back("2#mu2e ggF-like");
  categories->push_back("VBF"); color.push_back(kRed); marker.push_back(32); catlabels.push_back("VBF-like");
  cuts.push_back("(pass_vtx4lCut==1 && 140<m4l_constrained_HM&&m4l_constrained_HM<3000&&event_type==0 && !(dijet_invmass>400 && dijet_deltaeta>3.3))");
  cuts.push_back("(pass_vtx4lCut==1 && 140<m4l_constrained_HM&&m4l_constrained_HM<3000&&event_type==1 && !(dijet_invmass>400 && dijet_deltaeta>3.3))");
  cuts.push_back("(pass_vtx4lCut==1 && 140<m4l_constrained_HM&&m4l_constrained_HM<3000&&(event_type==3||event_type==2) && !(dijet_invmass>400 && dijet_deltaeta>3.3))");
  cuts.push_back("(pass_vtx4lCut==1 && 140<m4l_constrained_HM&&m4l_constrained_HM<3000 && (dijet_invmass>400 && dijet_deltaeta>3.3))");



  //LWA
  categories->push_back("ggF_4mu"); color.push_back(kViolet); marker.push_back(24); catlabels.push_back("4#mu ggF-like");
  categories->push_back("ggF_4e"); color.push_back(kOrange); marker.push_back(25); catlabels.push_back("4e ggF-like");
  categories->push_back("ggF_2mu2e"); color.push_back(kGreen+2); marker.push_back(26); catlabels.push_back("2#mu2e ggF-like");
  cuts.push_back("(pass_vtx4lCut==1 && 140<m4l_constrained_HM&&m4l_constrained_HM<1300&&event_type==0)");
  cuts.push_back("(pass_vtx4lCut==1 && 140<m4l_constrained_HM&&m4l_constrained_HM<1300&&event_type==1)");
  cuts.push_back("(pass_vtx4lCut==1 && 140<m4l_constrained_HM&&m4l_constrained_HM<1300&&(event_type==3||event_type==2))");

  //MELA
  categories->push_back("ggF_4mu_MELA1"); color.push_back(kMagenta+4); marker.push_back(20); catlabels.push_back("4#mu MELA 0.0-0.2");
  categories->push_back("ggF_4mu_MELA2"); color.push_back(kMagenta+3); marker.push_back(21); catlabels.push_back("4#mu MELA 0.2-0.4");
  categories->push_back("ggF_4mu_MELA3"); color.push_back(kMagenta+2); marker.push_back(22); catlabels.push_back("4#mu MELA 0.4-0.6");
  categories->push_back("ggF_4mu_MELA4"); color.push_back(kMagenta+1); marker.push_back(23); catlabels.push_back("4#mu MELA 0.6-0.8");
  categories->push_back("ggF_4mu_MELA5"); color.push_back(kMagenta); marker.push_back(24); catlabels.push_back("4#mu MELA 0.8-1.0");
  categories->push_back("ggF_4e_MELA1"); color.push_back(kCyan+4); marker.push_back(25); catlabels.push_back("4e MELA 0.0-0.2");
  categories->push_back("ggF_4e_MELA2"); color.push_back(kCyan+3); marker.push_back(26); catlabels.push_back("4e MELA 0.2-0.4");
  categories->push_back("ggF_4e_MELA3"); color.push_back(kCyan+2); marker.push_back(27); catlabels.push_back("4e MELA 0.4-0.6");
  categories->push_back("ggF_4e_MELA4"); color.push_back(kCyan+1); marker.push_back(28); catlabels.push_back("4e MELA 0.6-0.8");
  categories->push_back("ggF_4e_MELA5"); color.push_back(kCyan); marker.push_back(29); catlabels.push_back("4e MELA 0.8-1.0");
  categories->push_back("ggF_2mu2e_MELA1"); color.push_back(kRed+4); marker.push_back(30); catlabels.push_back("2mu2e MELA 0.0-0.2");
  categories->push_back("ggF_2mu2e_MELA2"); color.push_back(kRed+3); marker.push_back(31); catlabels.push_back("2mu2e MELA 0.2-0.4");
  categories->push_back("ggF_2mu2e_MELA3"); color.push_back(kRed+2); marker.push_back(32); catlabels.push_back("2mu2e MELA 0.4-0.6");
  categories->push_back("ggF_2mu2e_MELA4"); color.push_back(kRed+1); marker.push_back(33); catlabels.push_back("2mu2e MELA 0.6-0.8");
  categories->push_back("ggF_2mu2e_MELA5"); color.push_back(kRed); marker.push_back(34); catlabels.push_back("2mu2e MELA 0.8-1.0");
  //const char* MELA = "(1/(1 + c_factor(m4l_constrained_HM,event_type)*exp(-1*KD_discriminant)))";
  const char* MELA = "MELA";
  cuts.push_back(Form("((pass_vtx4lCut==1 && 200<m4l_constrained_HM&&m4l_constrained_HM<1200&&event_type==0) && %s>=0.0 && %s<0.2)",MELA,MELA));
  cuts.push_back(Form("((pass_vtx4lCut==1 && 200<m4l_constrained_HM&&m4l_constrained_HM<1200&&event_type==0) && %s>=0.2 && %s<0.4)",MELA,MELA));
  cuts.push_back(Form("((pass_vtx4lCut==1 && 200<m4l_constrained_HM&&m4l_constrained_HM<1200&&event_type==0) && %s>=0.4 && %s<0.6)",MELA,MELA));
  cuts.push_back(Form("((pass_vtx4lCut==1 && 200<m4l_constrained_HM&&m4l_constrained_HM<1200&&event_type==0) && %s>=0.6 && %s<0.8)",MELA,MELA));
  cuts.push_back(Form("((pass_vtx4lCut==1 && 200<m4l_constrained_HM&&m4l_constrained_HM<1200&&event_type==0) && %s>=0.8 && %s<=1.0)",MELA,MELA));
  cuts.push_back(Form("((pass_vtx4lCut==1 && 200<m4l_constrained_HM&&m4l_constrained_HM<1200&&event_type==1) && %s>=0.0 && %s<0.2)",MELA,MELA));
  cuts.push_back(Form("((pass_vtx4lCut==1 && 200<m4l_constrained_HM&&m4l_constrained_HM<1200&&event_type==1) && %s>=0.2 && %s<0.4)",MELA,MELA));
  cuts.push_back(Form("((pass_vtx4lCut==1 && 200<m4l_constrained_HM&&m4l_constrained_HM<1200&&event_type==1) && %s>=0.4 && %s<0.6)",MELA,MELA));
  cuts.push_back(Form("((pass_vtx4lCut==1 && 200<m4l_constrained_HM&&m4l_constrained_HM<1200&&event_type==1) && %s>=0.6 && %s<0.8)",MELA,MELA));
  cuts.push_back(Form("((pass_vtx4lCut==1 && 200<m4l_constrained_HM&&m4l_constrained_HM<1200&&event_type==1) && %s>=0.8 && %s<=1.0)",MELA,MELA));
  cuts.push_back(Form("((pass_vtx4lCut==1 && 200<m4l_constrained_HM&&m4l_constrained_HM<1200&&(event_type==2||event_type==3)) && %s>=0.0 && %s<0.2)",MELA,MELA));
  cuts.push_back(Form("((pass_vtx4lCut==1 && 200<m4l_constrained_HM&&m4l_constrained_HM<1200&&(event_type==2||event_type==3)) && %s>=0.2 && %s<0.4)",MELA,MELA));
  cuts.push_back(Form("((pass_vtx4lCut==1 && 200<m4l_constrained_HM&&m4l_constrained_HM<1200&&(event_type==2||event_type==3)) && %s>=0.4 && %s<0.6)",MELA,MELA));
  cuts.push_back(Form("((pass_vtx4lCut==1 && 200<m4l_constrained_HM&&m4l_constrained_HM<1200&&(event_type==2||event_type==3)) && %s>=0.6 && %s<0.8)",MELA,MELA));
  cuts.push_back(Form("((pass_vtx4lCut==1 && 200<m4l_constrained_HM&&m4l_constrained_HM<1200&&(event_type==2||event_type==3)) && %s>=0.8 && %s<=1.0)",MELA,MELA));
  */
