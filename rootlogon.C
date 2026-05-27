{
    cout << "=== Custom ROOT Environment Loading ===" << endl;
    cout << "Compiling PolHistograms.h..." << endl;
    gROOT->ProcessLine(".L PolHistograms.h+");
    
    cout << "Compiling PolAnalysis.h..." << endl;
    gROOT->ProcessLine(".L PolAnalysis.h+");
    
    cout << "Classes compiled and ready!" << endl;
    gROOT->ProcessLine(".L analyze_blpmc.C");
	cout << "Usage: .x analyze_blpmc.C" << endl;
    
    cout << "====================================" << endl;
}