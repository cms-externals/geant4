//
// Created by Stewart Boogert on 28/02/2023.
//

#ifndef G4VTKGRIDPIPELINE_HH
#define G4VTKGRIDPIPELINE_HH

class G4VtkStructuredGridPipeline
{
    G4VtkStructuredGridPipeline(G4int nxIn, G4int nyIn, G4int nzIn) : nx(nxIn), ny(nyIn), nz(nzIn)
    {
      structuredGrid->SetDimensions(nx, ny, nz);
      structuredGrid->SetPoints(points);
      structuredGrid->GetCellData()->SetScalars(cellValues);
      structuredGrid->GetPointData()->SetScalars(pointValues);

      mapper->SetInputData(structuredGrid);
      actor->SetMapper(mapper);
    }

    void Modified(){};
    void Clear(){};
    void Print(){};

    G4int nx, ny, nz;
    vtkSmartPointer<vtkPoints> points;
    vtkSmartPointer<vtkDoubleArray> pointValues;
    vtkSmartPointer<vtkDoubleArray> cellValues;
    vtkSmartPointer<vtkStructuredGrid> structuredGrid;
    vtkSmartPointer<vtkDataSetMapper> mapper;
    vtkSmartPointer<vtkActor> actor;
};

#endif  // G4VTKGRIDPIPELINE_HH
