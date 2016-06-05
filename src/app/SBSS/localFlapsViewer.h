/////////////////////////////////////////////////////////////////////////////
// Name:        localFlapsViewer.h
// Purpose:     wxWidgets executive for localFlaps program
// Author:      Court Cutting, MD
// Date:		March 11, 2012
// Copyright:   (c) Court Cutting - All rights reserved.
/////////////////////////////////////////////////////////////////////////////

#ifndef __LOCAL_FLAPS_VIEWER_H__
#define __LOCAL_FLAPS_VIEWER_H__

#include "wxGraphics.h"
#include <wx/glcanvas.h>
#include <wx/cmdline.h>
#include "surgicalActions.h"

// forward declarations
class flapToolsDlg;

// The OpenGL-enabled canvas
class BaseGLCanvas : public wxGLCanvas
{
public:
    BaseGLCanvas(wxWindow *parent,
                 wxWindowID id = wxID_ANY,
                 int *gl_attrib = NULL);

//    BaseGLCanvas(wxWindow *parent, wxWindowID id = wxID_ANY,
//        const wxPoint& pos = wxDefaultPosition,
//        const wxSize& size = wxDefaultSize, long style = 0,
//        const wxString& name = wxT("TestGLCanvas"),
//        int *attributes = 0);

    virtual ~BaseGLCanvas();

	void setSurgicalActions(surgicalActions *sa)	{_surgAct = sa;}
	wxGraphics* getWxGraphics() {return &_wGG;}

    void OnPaint(wxPaintEvent& event);
    void OnSize(wxSizeEvent& event);
    void OnKeyDown(wxKeyEvent& event);
    void OnKeyUp(wxKeyEvent& event);
    void OnMouseEvent(wxMouseEvent& event);
    void OnIdle(wxIdleEvent &event);
//	bool HandleKeyDown(WXWPARAM wParam, WXLPARAM lParam);


private:
	surgicalActions *_surgAct;
    wxGLContext* m_glRC;
	wxGraphics _wGG;
	MainFrame *_frame;
	bool _surgicalDrag,_glInitialized;
	int _lastSurgX,_lastSurgY,_xSize,_ySize;

    wxDECLARE_NO_COPY_CLASS(BaseGLCanvas);
    DECLARE_EVENT_TABLE()
};


// The frame containing the GL canvas
class MainFrame : public wxFrame
{
public:
    MainFrame(wxFrame *frame,
            const wxString& title,
            const wxPoint& pos = wxDefaultPosition,
            const wxSize& size = wxDefaultSize,
//            long style = wxDEFAULT_FRAME_STYLE|wxWANTS_CHARS);
            long style = wxMINIMIZE_BOX|wxMAXIMIZE_BOX|wxRESIZE_BORDER|wxSYSTEM_MENU|wxCLOSE_BOX);

    virtual ~MainFrame();

	void setToolState(int toolState);
	void setSurgicalActions(surgicalActions *surgAct) {_surgAct=surgAct; } //_pms = _surgAct->getPhysbamMayaScene();}
	surgicalActions* getSurgicalActions(){return _surgAct;}
	void OnToolsHook(wxCommandEvent& WXUNUSED(event) );
	void OnToolsScalpel(wxCommandEvent& WXUNUSED(event) );
	void OnToolsSuture(wxCommandEvent& WXUNUSED(event) );
	void OnToolsExcise(wxCommandEvent& WXUNUSED(event) );
	void OnToolsView(wxCommandEvent& WXUNUSED(event) );
	void OnModelViewWire(wxCommandEvent& WXUNUSED(event) );
	void OnModelViewTex(wxCommandEvent& WXUNUSED(event) );
	void OnModelViewNorm(wxCommandEvent& WXUNUSED(event) );
	void OnModelViewClear(wxCommandEvent& WXUNUSED(event) );
	void OnModelSmoothTexNorms(wxCommandEvent& WXUNUSED(event) );
	void OnModelScale(wxCommandEvent& WXUNUSED(event) );
	void OnModelCreateDirichletFile(wxCommandEvent& WXUNUSED(event) );
	void OnModelSaveChanged(wxCommandEvent& WXUNUSED(event) );
	void OnServerOff(wxCommandEvent& WXUNUSED(event) );
	void OnServerPassive(wxCommandEvent& WXUNUSED(event) );
	void OnServerActive(wxCommandEvent& WXUNUSED(event) );
	void toggleShowToolbox(wxCommandEvent& WXUNUSED(event) );
	void LoadModuleFile(wxCommandEvent& WXUNUSED(event) );
	void LoadHistoryFile(wxCommandEvent& WXUNUSED(event) );
	void nextHistoryEvent(wxCommandEvent& WXUNUSED(event) );
	void SaveHistoryFile(wxCommandEvent& WXUNUSED(event) );
	void getFilePath(const char *dialogTitle, const char *startPath, const char *fileFilter, bool openNotSave, std::string &returnDirectory, std::string &returnFilename);
	int showModalMessageDialog(const char *message, const char *messageTitle, bool OKnotYesNoButtons, bool cancelButton);
#if !defined(NO_PHYSICS) && defined(USE_RPC_INTERFACE)
    int getRPCServerDialog( std::string& server_address, std::string& server_port);
#endif
	BaseGLCanvas* getBaseGLCanvas() {return m_canvas;}
	bool IsCtrlOrShiftKeyDown();

    flapToolsDlg	*_ftd;
	surgicalActions *_surgAct;
    BaseGLCanvas *m_canvas;

private :
    void OnExit(wxCommandEvent& event);
    wxMenu *_toolsMenu;
    wxMenu *_serverMenu;
	int _toolState;
	bool _showToolbox;

    DECLARE_EVENT_TABLE()
};

// Define a new application type
class MyApp : public wxApp
{
public:
    virtual bool OnInit();
//	virtual int FilterEvent(wxEvent& event);
private:
	surgicalActions _surgAct;
	MainFrame *_frame;
	bool _cmdLineInput;
};

static const wxCmdLineEntryDesc g_cmdLineDesc [] =
{
     { wxCMD_LINE_OPTION , "f", "file", "scene file option",
          wxCMD_LINE_VAL_STRING , wxCMD_LINE_PARAM_OPTIONAL  },
     { wxCMD_LINE_OPTION , "h", "history", "history file option",
          wxCMD_LINE_VAL_STRING , wxCMD_LINE_PARAM_OPTIONAL  },
     { wxCMD_LINE_SWITCH, "a", "active", "active client",
          wxCMD_LINE_VAL_NONE, wxCMD_LINE_PARAM_OPTIONAL  },
     { wxCMD_LINE_SWITCH, "n", "noserver", "disables the server" },
     { wxCMD_LINE_SWITCH, "p", "passive", "passive client" },
 
     { wxCMD_LINE_NONE }
};

// Menu IDs
enum
{
    TOOLS_VIEW,
    TOOLS_HOOK,
    TOOLS_SCALPEL,
    TOOLS_SUTURE,
    TOOLS_EXCISE,
    TOOLS_SHOW_TOOLBOX,
    MODEL_VIEW_WIREFRAME,
    MODEL_VIEW_TEX_SEAMS,
    MODEL_VIEW_NORM_SEAMS,
    MODEL_VIEW_CLEAR,
    MODEL_SMOOTH_TEXTURES_NORMS,
    MODEL_SCALE,
    MODEL_DIRICHLET,
    SERVER_OFF,
    SERVER_PASSIVE,
    SERVER_ACTIVE,
    HISTORY_SAVE,
    HISTORY_LOAD,
    HISTORY_NEXT,
    MODEL_SAVE,
    MODULE_FILE_OPEN
};

#endif // __LOCAL_FLAPS_VIEWER_H__
