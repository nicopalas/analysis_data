void SaveCutG(char *name)
{
  TCutG *mycutg;
  mycutg = (TCutG *)gROOT->GetListOfSpecials()->FindObject("CUTG");
  mycutg->SetName(name);

  cout << "--------------------    SetCutG()   --------------------------" << endl;
  cout << " You have created a graphical cut with the following features:" << endl;
  cout << "  Name: " << mycutg->GetName() << endl;
  cout << "  Title: " << mycutg->GetTitle() << endl;
  cout << "  X variable: " << mycutg->GetVarX() << endl;
  cout << "  Y variable: " << mycutg->GetVarY() << endl;
  cout << "-------------------------------------------------------------" << endl;
}

void WriteCutG(char *name)
{
  TCutG *mycutg;
  mycutg = (TCutG *)gROOT->GetListOfSpecials()->FindObject(name);
  mycutg->Write();
}

void ReadCutG(char *name, char *nameOfFile)
{
  const char pwDirectory[50] = gDirectory->GetName();
  sprintf(pwDirectory, "%s%s", pwDirectory, ":/");

  TFile *thefile = new TFile(nameOfFile);
  TCutG *mycutg;
  thefile->GetObject(name, mycutg);
  gDirectory->cd(pwDirectory);

  cout << "------------------    ReadCutG()   ------------------------" << endl;
  cout << " You have read a graphical cut from the file " << nameOfFile
       << " with the following features:" << endl;
  cout << "  Name: " << mycutg->GetName() << endl;
  cout << "  Title: " << mycutg->GetTitle() << endl;
  cout << "  X variable: " << mycutg->GetVarX() << endl;
  cout << "  Y variable: " << mycutg->GetVarY() << endl;
  cout << "-------------------------------------------------------------" << endl;
}

void help()
{
  cout << "--------------------    help()   --------------------------" << endl;
  cout << "For a help creating graphical cuts, please type: helpSaveCutG()" << endl;
  cout << "For a help writing in a file a graphical cut, please type: helpWriteCutG()" << endl;
  cout << "For a help reading from a file an existing graphical cut, please type: helpReadCutG()" << endl;
  cout << "Other helpers are under construction   " << endl;
}

void helpSaveCutG()
{
  cout << "----------------    helpSaveCutG()   ----------------------" << endl;
  cout << " How to make, store and use a graphical cut:               " << endl;
  cout << " 1- DRAW the 2D histogram where you want to make the cut   " << endl;
  cout << "       example:   h101->Draw(\"Ntfi:Sst1i\") " << endl;
  cout << " 2- ADD the toolbal to your Canvas (if not present), by    " << endl;
  cout << "    clicking on    View  -> Toolbar                        " << endl;
  cout << " 3- CLICK on the last icon of the toolbar (scissors)       " << endl;
  cout << " 4- SELECT the points of the graphical cut, by a single    " << endl;
  cout << "    click in each vertex, and a double click in the last   " << endl;
  cout << "    vertex. The double click will close a polygon          " << endl;
  cout << " 5- CALL the funcion SaveCutG(\"name\")  where name is     " << endl;
  cout << "    a meaningful name you want to give to the graphical cut" << endl;
  cout << "       example:    SaveCutG(\"mySuperCut\");               " << endl;
  cout << " 6- USE the graphical cut in other histograms              " << endl;
  cout << "       example:    h101->Draw(\"Eventno\",\"supercut\")    " << endl;
  cout << "-----------------------------------------------------------" << endl;
}

void helpWriteCutG()
{
  cout << "----------------    helpWriteCutG()   ---------------------" << endl;
  cout << " How to write your graphical cut in a root file            " << endl;
  cout << " 1 - Create your graphical cut as explained before         " << endl;
  cout << "     (type helpSaveCutG() in case of doubt)                " << endl;
  cout << " 2 - Create or open a file for storing your graphical cut: " << endl;
  cout << "   example: TFile myF(\"myCuts.root\",\"CREATE\",\"myCuts.root\")" << endl;
  cout << " 3 - Call the function (that is, type):   WriteCutG(\"namec\");" << endl;
  cout << "     where \"namec\" is the name of the graphical cut       " << endl;
  cout << "-----------------------------------------------------------" << endl;
}

void helpReadCutG()
{
  cout << "----------------    helpReadCutG()   ---------------------" << endl;
  cout << " How to read your graphical cut from a root file ...      " << endl;
  cout << " ...supposed you've created your graphical cut and stored " << endl;
  cout << " it as shown in  helpSaveCutG()  and  helpWriteCutG()     " << endl;
  cout << " 1 - Open the tree where you want to apply the cuts       " << endl;
  cout << " 2 - Call the function (that is, type): ReadCutG(namec,namef);" << endl;
  cout << "     where \"namec\" is the name of the graphical cut       " << endl;
  cout << "     and where \"namef\" is the name of the file where the cutg is" << endl;
  cout << " 3 - Now you can apply the cuts...                         " << endl;
  cout << "           example: h101->Draw(\"Eventno\",\"acut\")       " << endl;
  cout << "-----------------------------------------------------------" << endl;
}