/////////////////////////////////////////////////////////////////////////////
// Name:        flaptoolsdlg.cpp
// Purpose:     
// Author:      Court Cutting, M.D.
// Modified by: 
// Created:     02/10/2007 17:21:27
// RCS-ID:      
// Copyright:   All copyright rights reserved.
// Licence:     
/////////////////////////////////////////////////////////////////////////////

// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifdef __BORLANDC__
#pragma hdrstop
#endif

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

////@begin includes
#include <wx/event.h>
#include "wx/statline.h"
#include "surgicalActions.h"
#include "localFlapsViewer.h"
#include "flaptoolsdlg.h"
////@end includes


////@begin XPM images
////@end XPM images


/*!
 * flapToolsDlg type definition
 */

IMPLEMENT_DYNAMIC_CLASS( flapToolsDlg, wxDialog )


/*!
 * flapToolsDlg event table definition
 */

/* BEGIN_EVENT_TABLE( flapToolsDlg, wxDialog )
////@begin flapToolsDlg event table entries
    EVT_INIT_DIALOG( flapToolsDlg::OnInitDialog )
	EVT_RADIOBOX(ID_RADIOBOX_Tools,flapToolsDlg::OnToolChange )
	EVT_BUTTON(ID_BUTTON_Next,flapToolsDlg::OnHistoryNext )

////@end flapToolsDlg event table entries
END_EVENT_TABLE() */


/*!
 * flapToolsDlg constructors
 */

flapToolsDlg::flapToolsDlg()
{
    Init();
}

flapToolsDlg::flapToolsDlg( wxWindow* parent, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style )
{
    Init();
    Create(parent, id, caption, pos, size, style);
    Bind(wxEVT_INIT_DIALOG, &flapToolsDlg::OnInitDialog, this);
    Bind(wxEVT_COMMAND_RADIOBOX_SELECTED, &flapToolsDlg::OnToolChange, this, ID_RADIOBOX_Tools);
    Bind(wxEVT_COMMAND_BUTTON_CLICKED, &flapToolsDlg::OnHistoryNext, this, ID_BUTTON_Next);
    Bind(wxEVT_LEAVE_WINDOW, &flapToolsDlg::OnMouseLeavesWindow, this);
    Bind(wxEVT_CLOSE_WINDOW, &flapToolsDlg::OnCloseToolbox, this);
//    Bind(wxEVT(wxEVT_CLOSE_WINDOW, &flapToolsDlg::OnCloseToolbox, this);
//	void flapToolsDlg::OnKeyDown(wxKeyEvent& event)


}


/*!
 * flapToolsDlg creator
 */

bool flapToolsDlg::Create( wxWindow* parent, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style )
{
////@begin flapToolsDlg creation
    SetExtraStyle(wxWS_EX_BLOCK_EVENTS);
	wxPoint newPos(pos);
	newPos.y += 20;
//    wxDialog::Create( parent, id, caption, newPos, size, style );
    wxDialog::Create( parent, id, wxString(), newPos, size, style );

    CreateControls();
    if (GetSizer())
    {
        GetSizer()->SetSizeHints(this);
    }
    Centre();
////@end flapToolsDlg creation
	_parentFrame = (MainFrame*)parent;
    return true;
}


/*!
 * flapToolsDlg destructor
 */

flapToolsDlg::~flapToolsDlg()
{
////@begin flapToolsDlg destruction
////@end flapToolsDlg destruction
}


/*!
 * Member initialisation
 */

void flapToolsDlg::Init()
{
////@begin flapToolsDlg member initialisation
////@end flapToolsDlg member initialisation
}


/*
 * Control creation for flapToolsDlg
 */

void flapToolsDlg::CreateControls()
{    
////@begin flapToolsDlg content construction
    flapToolsDlg* itemDialog1 = this;
    wxBoxSizer* itemBoxSizer2 = new wxBoxSizer(wxVERTICAL);
    itemDialog1->SetSizer(itemBoxSizer2);
    wxArrayString itemRadioBox3Strings;
    itemRadioBox3Strings.Add(_("&View"));
    itemRadioBox3Strings.Add(_("&Hook"));
    itemRadioBox3Strings.Add(_("&Knife"));
    itemRadioBox3Strings.Add(_("&Suture"));
    itemRadioBox3Strings.Add(_("&Excise"));
    _itemRadioBox3 = new wxRadioBox( itemDialog1, ID_RADIOBOX_Tools, _("TOOLS"), wxDefaultPosition, wxDefaultSize, itemRadioBox3Strings, 5, wxRA_SPECIFY_ROWS );
    _itemRadioBox3->SetSelection(0);
    itemBoxSizer2->Add(_itemRadioBox3, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    wxStaticLine* itemStaticLine4 = new wxStaticLine( itemDialog1, ID_STATICLINE, wxDefaultPosition, wxDefaultSize, wxLI_HORIZONTAL );
    itemBoxSizer2->Add(itemStaticLine4, 0, wxGROW|wxALL, 5);

    wxButton* itemButton5 = new wxButton( itemDialog1, ID_BUTTON_Next, _("NEXT"), wxDefaultPosition, wxDefaultSize, 0 );
    itemBoxSizer2->Add(itemButton5, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);
////@end flapToolsDlg content construction
}


/*!
 * Should we show tooltips?
 */

bool flapToolsDlg::ShowToolTips()
{
    return true;
}

/*!
 * Get bitmap resources
 */

wxBitmap flapToolsDlg::GetBitmapResource( const wxString& name )
{
    // Bitmap retrieval
////@begin flapToolsDlg bitmap retrieval
    wxUnusedVar(name);
    return wxNullBitmap;
////@end flapToolsDlg bitmap retrieval
}

/*!
 * Get icon resources
 */

wxIcon flapToolsDlg::GetIconResource( const wxString& name )
{
    // Icon retrieval
////@begin flapToolsDlg icon retrieval
    wxUnusedVar(name);
    return wxNullIcon;
////@end flapToolsDlg icon retrieval
}


/*!
 * wxEVT_INIT_DIALOG event handler for ID_FLAPTOOLSDLG
 */

void flapToolsDlg::OnInitDialog( wxInitDialogEvent& event )
{
////@begin wxEVT_INIT_DIALOG event handler for ID_FLAPTOOLSDLG in flapToolsDlg.
    // Before editing this code, remove the block markers.
    event.Skip();
////@end wxEVT_INIT_DIALOG event handler for ID_FLAPTOOLSDLG in flapToolsDlg. 
}

void flapToolsDlg::OnToolChange( wxCommandEvent& event)
{
	// COURT - next line needs to be there or text control won't return to parent window
    event.Skip();
	int radioButtonNumber=_itemRadioBox3->GetSelection();
	switch (radioButtonNumber)
	{
	case  0 :
		_parentFrame->_surgAct->setToolState(0);
		_parentFrame->setToolState(0);
		break;
	case  1 :
		_parentFrame->_surgAct->setToolState(1);
		_parentFrame->setToolState(1);
		break;
	case  2 :
		_parentFrame->_surgAct->setToolState(2);
		_parentFrame->setToolState(2);
		break;
	case  3 :
		_parentFrame->_surgAct->setToolState(3);
		_parentFrame->setToolState(3);
		break;
	case  4 :
		_parentFrame->_surgAct->setToolState(4);
		_parentFrame->setToolState(4);
		break;
	default :
		_parentFrame->_surgAct->setToolState(0);
		_parentFrame->setToolState(0);
		break;
	}
	_parentFrame->Refresh(false);
}

void flapToolsDlg::setToolBoxSelection(int toolState)
{
	 _itemRadioBox3->SetSelection(toolState);
}

void flapToolsDlg::OnHistoryNext( wxCommandEvent& WXUNUSED(event))
{
	_parentFrame->_surgAct->nextHistoryAction();
	_parentFrame->Refresh(false);
}

void flapToolsDlg::OnCloseToolbox( wxCloseEvent& WXUNUSED(event))
{
	wxCommandEvent junkEvent;
	_parentFrame->toggleShowToolbox(junkEvent);
}

void flapToolsDlg::OnMouseLeavesWindow( wxMouseEvent& WXUNUSED(event))
{
	_parentFrame->getBaseGLCanvas()->SetFocus();
}

void flapToolsDlg::OnKeyDown(wxKeyEvent& event)
{
	_parentFrame->getBaseGLCanvas()->SetFocus();
}
