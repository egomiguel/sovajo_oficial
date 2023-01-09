#ifndef VTK_RESLICE_CALLBACK_H
#define VTK_RESLICE_CALLBACK_H

#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "vtkCommand.h"
#include "vtkResliceCursorWidget.h"
#include "vtkResliceCursorLineRepresentation.h"
#include "vtkPlaneSource.h"
#include "vtkImagePlaneWidget.h"


class vtkResliceCursorCallback : public vtkCommand
{
public:
	static vtkResliceCursorCallback *New()
	{
		return new vtkResliceCursorCallback;
	}
	void Execute(vtkObject *caller, unsigned long, void *callData) override;
	vtkResliceCursorCallback() {}
	vtkImagePlaneWidget* IPW[3];
	vtkResliceCursorWidget *RCW[3];
};

#endif