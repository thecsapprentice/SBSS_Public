/////////////////////////////////////////////////////////////////////////////
// Name:        flaptoolsdlg.h
// Purpose:     
// Author:      Court Cutting, M.D.
// Modified by: 
// Created:     02/10/2007 17:21:27
// RCS-ID:      
// Copyright:   All copyright rights reserved.
// License:
/////////////////////////////////////////////////////////////////////////////

#ifndef _FLAPTOOLSDLG_H_
#define _FLAPTOOLSDLG_H_

/*!
 * Includes
 */

////@begin includes
////@end includes

/*!
 * Forward declarations
 */

////@begin forward declarations
class MainFrame;
////@end forward declarations

/*!
 * Control identifiers
 */

////@begin control identifiers
#define ID_FLAPTOOLSDLG 10000
#define ID_RADIOBOX_Tools 10004
#define ID_STATICLINE 10006
#define ID_BUTTON_Next 10005

#ifdef WIN32
#define SYMBOL_FLAPTOOLSDLG_STYLE wxDEFAULT_DIALOG_STYLE|wxWANTS_CHARS|wxFRAME_FLOAT_ON_PARENT|wxSTAY_ON_TOP	// wxSTAY_ON_TOP for MS_Windows
#else
//#define SYMBOL_FLAPTOOLSDLG_STYLE wxCAPTION|wxRESIZE_BORDER|wxFRAME_FLOAT_ON_PARENT|wxNO_BORDER
#define SYMBOL_FLAPTOOLSDLG_STYLE wxCAPTION|wxRESIZE_BORDER|wxFRAME_FLOAT_ON_PARENT|wxCLOSE_BOX
#endif

//#define SYMBOL_FLAPTOOLSDLG_STYLE wxDEFAULT_DIALOG_STYLE|wxWANTS_CHARS|wxFRAME_FLOAT_ON_PARENT
#define SYMBOL_FLAPTOOLSDLG_TITLE _("TOOLS")
#define SYMBOL_FLAPTOOLSDLG_IDNAME ID_FLAPTOOLSDLG
#define SYMBOL_FLAPTOOLSDLG_SIZE wxDefaultSize			// wxSize(300, 300)
#define SYMBOL_FLAPTOOLSDLG_POSITION wxDefaultPosition
////@end control identifiers


/*!
 * flapToolsDlg class declaration
 */

class flapToolsDlg: public wxDialog
{    
    DECLARE_DYNAMIC_CLASS( flapToolsDlg )
//    DECLARE_EVENT_TABLE()

public:
    /// Constructors
    flapToolsDlg();
    flapToolsDlg( wxWindow* parent, wxWindowID id = SYMBOL_FLAPTOOLSDLG_IDNAME, const wxString& caption = SYMBOL_FLAPTOOLSDLG_TITLE, const wxPoint& pos = SYMBOL_FLAPTOOLSDLG_POSITION, const wxSize& size = SYMBOL_FLAPTOOLSDLG_SIZE, long style = SYMBOL_FLAPTOOLSDLG_STYLE );

    /// Creation
    bool Create( wxWindow* parent, wxWindowID id = SYMBOL_FLAPTOOLSDLG_IDNAME, const wxString& caption = SYMBOL_FLAPTOOLSDLG_TITLE, const wxPoint& pos = SYMBOL_FLAPTOOLSDLG_POSITION, const wxSize& size = SYMBOL_FLAPTOOLSDLG_SIZE, long style = SYMBOL_FLAPTOOLSDLG_STYLE );

    /// Destructor
    ~flapToolsDlg();

    /// Initialises member variables
    void Init();

    /// Creates the controls and sizers
    void CreateControls();

////@begin flapToolsDlg event handler declarations

    /// wxEVT_INIT_DIALOG event handler for ID_FLAPTOOLSDLG
    void OnInitDialog( wxInitDialogEvent& event );
	void OnToolChange( wxCommandEvent& event);
	void OnHistoryNext( wxCommandEvent& WXUNUSED(event));
	void OnMouseLeavesWindow( wxMouseEvent& WXUNUSED(event));
	void OnCloseToolbox( wxCloseEvent& WXUNUSED(event));
	void OnKeyDown(wxKeyEvent& event);

////@end flapToolsDlg event handler declarations

////@begin flapToolsDlg member function declarations

    /// Retrieves bitmap resources
    wxBitmap GetBitmapResource( const wxString& name );

    /// Retrieves icon resources
    wxIcon GetIconResource( const wxString& name );
	void setToolBoxSelection(int toolState);
////@end flapToolsDlg member function declarations

    /// Should we show tooltips?
    static bool ShowToolTips();

////@begin flapToolsDlg member variables
private:
	MainFrame* _parentFrame;
    wxRadioBox* _itemRadioBox3;
////@end flapToolsDlg member variables
};

#endif
    // _FLAPTOOLSDLG_H_
