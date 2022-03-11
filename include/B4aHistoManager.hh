#ifndef B4aHistoManager_h
#define B4aHistoManager_h 1

#include "globals.hh"

#include "g4root.hh"
//#include "g4xml.hh"

const G4int MaxNtCol1=50;

class HistoManager
{
  public:
   HistoManager();
  ~HistoManager();



    void Book();
    void save();
    void FillNtuple(G4int, G4double);
    void AddRow(G4int);

private:
  G4String fFileName;
  G4int MaxNtCol;
  G4int fNtColId[MaxNtCol1];
};



#endif
