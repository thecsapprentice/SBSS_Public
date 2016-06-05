/////////////////////////////////////////////////////////////////////////////
// Name:        localFlapsViewer.h
// Purpose:     wxWidgets executive for localFlaps program
// Author:      Court Cutting, MD
// Date:		March 11, 2012
// Copyright:   (c) Court Cutting - All rights reserved.
/////////////////////////////////////////////////////////////////////////////

// For compilers that support precompilation, includes "wx.h".
#include "wx/wxprec.h"

#ifdef __BORLANDC__
#pragma hdrstop
#endif

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

#if !wxUSE_GLCANVAS
    #error "OpenGL required: set wxUSE_GLCANVAS to 1 and rebuild the library"
#endif

#include "wx/timer.h"
//#include "wx/glcanvas.h"
#include "wx/math.h"
#include "wx/log.h"
#include "wx/cmdline.h"
#include "wx/wfstream.h"
#include "wx/zstream.h"
#include "wx/txtstrm.h"
#include "wx/filedlg.h"
#include "wx/filename.h"
#include "wx/filefn.h"
#include "wx/numdlg.h"
#include "wx/menu.h"
#include <sstream>

#include "flaptoolsdlg.h"
#include "localFlapsViewer.h"

#if !defined(NO_PHYSICS) && defined(USE_RPC_INTERFACE)
#include <wx/uri.h>
#include <wx/textctrl.h>
#endif

#if !wxUSE_MENUS
    // nice try...
    #error "menu sample requires wxUSE_MENUS=1"
#endif // wxUSE_MENUS

// not all ports have support for EVT_CONTEXT_MENU yet, don't define
// USE_CONTEXT_MENU for those which don't
#if defined(__WXMOTIF__) || defined(__WXPM__) || defined(__WXX11__)
    #define USE_CONTEXT_MENU 0
#else
    #define USE_CONTEXT_MENU 1
#endif

// this sample is useful when a new port is developed
// and usually a new port has majority of flags turned off
#if wxUSE_LOG && wxUSE_TEXTCTRL
    #define USE_LOG_WINDOW 1
#else
    #define USE_LOG_WINDOW 0
#endif

#if !defined(NO_PHYSICS) && !defined(USE_RPC_INTERFACE)
#include <Thread_Queueing/PTHREAD_QUEUE.h>

using namespace PhysBAM;
extern PTHREAD_QUEUE* pthread_queue;

#endif



//---------------------------------------------------------------------------
// MyApp
//---------------------------------------------------------------------------

IMPLEMENT_APP(MyApp)

// `Main program' equivalent, creating windows and returning main app frame
bool MyApp::OnInit()
{
    // parse comand line
    wxString cmdSceneFile,cmdHistoryFile;
    wxCmdLineParser cmdParser(g_cmdLineDesc, argc, argv);
    int serverMode=-1,res,fileLoad=0;
    wxLogNull log;
    // false suppresses AutoUsage message
    res = cmdParser.Parse(false);
    if(res==-1 || res>0 || cmdParser.Found(wxT("h"))) {
        cmdParser.Usage();
        return(false);
    }
    if(cmdParser.Found(wxT("n")))
        serverMode = 0;
    if(cmdParser.Found(wxT("p")))
        serverMode = 1;
    if(cmdParser.Found(wxT("a")))
        serverMode = 2;
    if(cmdParser.Found(wxT("f"),&cmdSceneFile))
        fileLoad |= 1;
    if(cmdParser.Found(wxT("h"),&cmdHistoryFile))
        fileLoad |= 2;

//    if ( !wxApp::OnInit() )
//        return false;

    #if !defined(NO_PHYSICS) && !defined(USE_RPC_INTERFACE)

    pthread_queue=new PhysBAM::PTHREAD_QUEUE(12);  

    #endif


#if !defined(NO_PHYSICS) && defined(USE_RPC_INTERFACE)
    std::string server;
    std::string port;
    _frame->getRPCServerDialog(server, port);
    _surgAct._cle.InitializeCLERPC( server, port);
#endif


    // Create the main frame window
	_frame = new MainFrame(NULL, wxT("Local Flaps Window"),wxDefaultPosition, wxDefaultSize, wxDEFAULT_FRAME_STYLE);
	_frame->setSurgicalActions(&_surgAct);
	_frame->m_canvas->setSurgicalActions(&_surgAct);
	_surgAct.setWxGraphics(_frame->m_canvas->getWxGraphics());
	_surgAct.setMainFrame(_frame);
	// Show the frame
	_frame->Show(true);
	if(serverMode<1)
		_frame->Maximize();
	if(serverMode>-1)
		_surgAct.setServerMode(serverMode);
	if(fileLoad) {
		wxString inFile;
		if(fileLoad>1)
			inFile = cmdHistoryFile;
		else
			inFile = cmdSceneFile;
		wxFileName fName(inFile);
		fName.Normalize(wxPATH_NORM_LONG|wxPATH_NORM_DOTS|wxPATH_NORM_TILDE|wxPATH_NORM_ABSOLUTE);
		if(fileLoad>1)
			_surgAct.loadHistory(fName.GetFullPath().c_str());
		else
			_surgAct.loadScene(fName.GetPathWithSep().c_str(),fName.GetFullName().c_str(),true);
	}

    return true;
}

//---------------------------------------------------------------------------
// MyFrame
//---------------------------------------------------------------------------

BEGIN_EVENT_TABLE(MainFrame, wxFrame)
    EVT_MENU(wxID_EXIT, MainFrame::OnExit)
    EVT_MENU(MODULE_FILE_OPEN,MainFrame::LoadModuleFile)
    EVT_MENU(TOOLS_VIEW,MainFrame::OnToolsView)
    EVT_MENU(TOOLS_HOOK,MainFrame::OnToolsHook)
    EVT_MENU(TOOLS_SCALPEL,MainFrame::OnToolsScalpel)
    EVT_MENU(TOOLS_SUTURE,MainFrame::OnToolsSuture)
    EVT_MENU(TOOLS_EXCISE,MainFrame::OnToolsExcise)
    EVT_MENU(TOOLS_SHOW_TOOLBOX,MainFrame::toggleShowToolbox)
    EVT_MENU(HISTORY_SAVE,MainFrame::SaveHistoryFile)
    EVT_MENU(HISTORY_LOAD,MainFrame::LoadHistoryFile)
    EVT_MENU(HISTORY_NEXT,MainFrame::nextHistoryEvent)
    EVT_MENU(MODEL_VIEW_WIREFRAME,MainFrame::OnModelViewWire)
    EVT_MENU(MODEL_VIEW_TEX_SEAMS,MainFrame::OnModelViewTex)
    EVT_MENU(MODEL_VIEW_NORM_SEAMS,MainFrame::OnModelViewNorm)
    EVT_MENU(MODEL_VIEW_CLEAR,MainFrame::OnModelViewClear)
    EVT_MENU(MODEL_SMOOTH_TEXTURES_NORMS,MainFrame::OnModelSmoothTexNorms)
    EVT_MENU(MODEL_SCALE,MainFrame::OnModelScale)
    EVT_MENU(MODEL_DIRICHLET,MainFrame::OnModelCreateDirichletFile)
    EVT_MENU(MODEL_SAVE,MainFrame::OnModelSaveChanged)
    EVT_MENU(SERVER_OFF,MainFrame::OnServerOff)
    EVT_MENU(SERVER_PASSIVE,MainFrame::OnServerPassive)
    EVT_MENU(SERVER_ACTIVE,MainFrame::OnServerActive)
END_EVENT_TABLE()

void MainFrame::OnServerOff(wxCommandEvent& WXUNUSED(event) )
{
	_surgAct->setServerMode(0);
}

void MainFrame::OnServerPassive(wxCommandEvent& WXUNUSED(event) )
{
	_surgAct->setServerMode(1);
}

void MainFrame::OnServerActive(wxCommandEvent& WXUNUSED(event) )
{
	_surgAct->setServerMode(2);
}

void MainFrame::OnModelViewWire(wxCommandEvent& WXUNUSED(event) )
{
	_surgAct->showWireModel(3);
}

void MainFrame::OnModelViewTex(wxCommandEvent& WXUNUSED(event) )
{
	_surgAct->showWireModel(2);
}

void MainFrame::OnModelViewNorm(wxCommandEvent& WXUNUSED(event) )
{
	_surgAct->showWireModel(1);
}

void MainFrame::OnModelViewClear(wxCommandEvent& WXUNUSED(event) )
{
	_surgAct->showWireModel(0);
}

void MainFrame::OnModelSmoothTexNorms(wxCommandEvent& WXUNUSED(event) )
{
//	int tnn;
//	_surgAct->_cle.getIncisionPointer()->cleanModel(_surgAct->_cle.getElasticGlslTriangle(),0.005f,0.45f,tnn);
}

void MainFrame::OnModelScale(wxCommandEvent& WXUNUSED(event) )
{
	float maxDim=-1.0f,bbox[6];
	_surgAct->_cle.getBoundingBox(bbox);
	if(bbox[1]-bbox[0]>bbox[3]-bbox[2])
		maxDim = bbox[1]-bbox[0];
	else
		maxDim = bbox[3]-bbox[2];
	if(maxDim<bbox[5]-bbox[4])
		maxDim = bbox[5]-bbox[4];
	std::ostringstream istr;
	istr << "Current maximum cardinal axis dimension is " << maxDim << " meters.";
	wxString msg(istr.str().c_str(),wxConvUTF8);
	long value;
	wxNumberEntryDialog numDlg(this,msg,wxString("Enter a percent scale factor:"),wxString("Scale Factor Entry"),100L,1L,30000L);
    numDlg.CentreOnParent();
	if(numDlg.ShowModal()==wxID_OK)	{
		value = numDlg.GetValue();
		_surgAct->_cle.scaleScene(value*0.01f);
	}
}

void MainFrame::OnModelCreateDirichletFile(wxCommandEvent& WXUNUSED(event) )
{
	_surgAct->_cle.createFixedNodeFile();
}

void MainFrame::OnModelSaveChanged(wxCommandEvent& WXUNUSED(event) )
{
	_surgAct->_cle.saveModifiedScene();
}

bool MainFrame::IsCtrlOrShiftKeyDown()
{
	if (wxGetKeyState(WXK_CONTROL) || wxGetKeyState(WXK_SHIFT))
		return true;
	else
		return false;
}

void MainFrame::OnToolsHook(wxCommandEvent& WXUNUSED(event) )
{
	setToolState(1);
}

void MainFrame::OnToolsScalpel(wxCommandEvent& WXUNUSED(event) )
{
	setToolState(2);
}

void MainFrame::OnToolsSuture(wxCommandEvent& WXUNUSED(event) )
{
	setToolState(3);
}

void MainFrame::OnToolsExcise(wxCommandEvent& WXUNUSED(event) )
{
	setToolState(4);
}

void MainFrame::OnToolsView(wxCommandEvent& WXUNUSED(event) )
{
	setToolState(0);
}

void MainFrame::toggleShowToolbox(wxCommandEvent& WXUNUSED(event) )
{
	if(_showToolbox)	{
		_showToolbox=false;
		_toolsMenu->Check(TOOLS_SHOW_TOOLBOX,false);
		_ftd->Hide();
	}
	else	{
		_showToolbox=true;
		_toolsMenu->Check(TOOLS_SHOW_TOOLBOX,true);
		_ftd->Show();
	}
}

// My frame constructor
MainFrame::MainFrame(wxFrame *frame, const wxString& title, const wxPoint& pos,
    const wxSize& size, long style)
    : wxFrame(frame, wxID_ANY, title, pos, size, style),
      m_canvas(NULL)
{
//    SetIcon(wxICON(sample));
  // JACS
#ifdef __WXMSW__
    int *gl_attrib = NULL;
#else
    int gl_attrib[20] =
        { WX_GL_RGBA, WX_GL_MIN_RED, 1, WX_GL_MIN_GREEN, 1,
        WX_GL_MIN_BLUE, 1, WX_GL_DEPTH_SIZE, 1,
        WX_GL_DOUBLEBUFFER,
#  if defined(__WXMAC__) || defined(__WXCOCOA__)
        GL_NONE };
#  else
        None };
#  endif
#endif

	wxPoint framePos = GetPosition();
	framePos.x += 5;
	framePos.y += 55;
//	this->SetFocusIgnoringChildren();	// keyboard focus stays here and doesn't shift to dialog. Not a member of wx3.0.0
	_ftd = new flapToolsDlg(this);
	_ftd->SetPosition(framePos);
	_ftd->Show();
		// Make a menubar
    wxMenu *moduleMenu = new wxMenu;
    moduleMenu->Append(MODULE_FILE_OPEN,  wxT("&Load module"));
    moduleMenu->AppendSeparator();
    moduleMenu->Append(wxID_EXIT, wxT("E&xit"));
    moduleMenu->AppendSeparator();
    moduleMenu->Append(HISTORY_SAVE,  wxT("&Save History File"));
    moduleMenu->AppendSeparator();
    moduleMenu->AppendRadioItem(SERVER_OFF,  wxT("Server &Off"));
    moduleMenu->AppendRadioItem(SERVER_PASSIVE,  wxT("&Passive Client"));
    moduleMenu->AppendRadioItem(SERVER_ACTIVE,  wxT("&Active Client"));

    wxMenu *historyMenu = new wxMenu;
    historyMenu->Append(HISTORY_SAVE,  wxT("&Save"));
    historyMenu->Append(HISTORY_LOAD,  wxT("L&oad"));
    historyMenu->Append(HISTORY_NEXT,  wxT("&Next"));

    _toolsMenu = new wxMenu;	// this one is a saved member variable as we go back to check it frequently
    _toolsMenu->AppendRadioItem(TOOLS_VIEW,  wxT("&View"));
    _toolsMenu->AppendRadioItem(TOOLS_HOOK,  wxT("&Hook"));
    _toolsMenu->AppendRadioItem(TOOLS_SCALPEL,  wxT("&Scalpel"));
    _toolsMenu->AppendRadioItem(TOOLS_SUTURE,  wxT("S&uture"));
    _toolsMenu->AppendRadioItem(TOOLS_EXCISE,  wxT("&Excise"));
    _toolsMenu->AppendSeparator();
    _toolsMenu->AppendCheckItem(TOOLS_SHOW_TOOLBOX,  wxT("Show&Toolbox"));

    wxMenu *modelMenu = new wxMenu;
    modelMenu->Append(MODEL_VIEW_WIREFRAME,  wxT("view triangle wireframe"));
    modelMenu->Append(MODEL_VIEW_TEX_SEAMS,  wxT("view texture seam vertices"));
    modelMenu->Append(MODEL_VIEW_NORM_SEAMS,  wxT("view normal seam vertices"));
    modelMenu->Append(MODEL_VIEW_CLEAR,  wxT("clear wire view"));
    modelMenu->Append(MODEL_SMOOTH_TEXTURES_NORMS,  wxT("smooth texture and normal seam vertices"));
    modelMenu->Append(MODEL_SCALE,  wxT("scale model"));
    modelMenu->Append(MODEL_DIRICHLET,  wxT("make fixed geometry file"));
    modelMenu->Append(MODEL_SAVE,  wxT("save changed model"));

    _serverMenu = new wxMenu;
    _serverMenu->AppendRadioItem(SERVER_OFF,  wxT("&Off"));
    _serverMenu->AppendRadioItem(SERVER_PASSIVE,  wxT("&Passive client"));
    _serverMenu->AppendRadioItem(SERVER_ACTIVE,  wxT("&Active client"));

    wxMenuBar *menubar = new wxMenuBar( wxMB_DOCKABLE );
    bool success = menubar->Append(moduleMenu, wxT("&Module"));
    success = menubar->Append(historyMenu, wxT("&History"));
    success = menubar->Append(_toolsMenu, _T("&Tools"));
    success = menubar->Append(modelMenu, _T("&Model"));
    success = menubar->Append(_serverMenu, _T("&Server"));
    _showToolbox = true;

    SetMenuBar(menubar);
    _toolsMenu->Check(TOOLS_VIEW,true);
    _toolsMenu->Check(TOOLS_SHOW_TOOLBOX,true);
    _serverMenu->Check(SERVER_OFF,true);

	m_canvas = new BaseGLCanvas(this, wxID_ANY, gl_attrib);
	// next line has no effect. Toolbox on top seems to give it keyboard focus which is unshakeable.
//	m_canvas->SetFocus();	// sets keyboard focus on this window, not the tool dialog. Is a wxWindow call to base class
}

MainFrame::~MainFrame()
{
    delete m_canvas;
}

// Intercept menu commands
void MainFrame::OnExit( wxCommandEvent& WXUNUSED(event) )
{
//	_ftd->Destroy();	// COURT - I thought this was necessary, but it causes error "trying to close something that doesn't think is its parent"
    // true is to force the frame to close
    Close(true);
}

void MainFrame::LoadModuleFile(wxCommandEvent& WXUNUSED(event) )
{
	std::string retDir,retFname;
	getFilePath("Choose a surgical module file-","./","*.smd",true,retDir,retFname);
	if(retFname=="")
		return;
	_surgAct->loadScene(retDir.c_str(),retFname.c_str(),true);
    m_canvas->Refresh(false);

}

void MainFrame::LoadHistoryFile(wxCommandEvent& WXUNUSED(event) )
{
	if(!_surgAct->historyEmpty())	{
		showModalMessageDialog("Can only load a history file at the beginning. Some surgical history is already present.","USER HISTORY ERROR",true,false);
		return;
	}
	std::string retPath,retDir;
	getFilePath("Choose a surgical history file-","./Histories","*.hst",true,retDir,retPath);
	retDir += retPath;
	_surgAct->loadHistory(retDir.c_str());
}

void MainFrame::nextHistoryEvent(wxCommandEvent& WXUNUSED(event) )
{
	_surgAct->nextHistoryAction();
}

void MainFrame::SaveHistoryFile(wxCommandEvent& WXUNUSED(event) )
{
	std::string retPath,retDir;
	getFilePath("Name a history file which will save your surgical procedure-","./Histories","*.hst",false,retDir,retPath);
	size_t len = retPath.rfind(".hst");
	if(len>retPath.size())
		retPath.append(".hst");
	retDir += retPath;
	_surgAct->saveSurgicalHistory(retDir.c_str());
}


void MainFrame::getFilePath(const char *dialogTitle, const char *startPath, const char *fileFilter, bool openNotSave, std::string &returnDirectory, std::string &returnFilename)
{
	wxString s1(dialogTitle,wxConvUTF8);
	wxString s2(startPath,wxConvUTF8);
	wxString s3(fileFilter,wxConvUTF8);
	long styleFlags;
	if(openNotSave)
		styleFlags = wxFD_OPEN|wxFD_FILE_MUST_EXIST;
	else
		styleFlags = wxFD_SAVE|wxFD_OVERWRITE_PROMPT;
    wxFileDialog dialog
	(
		this,
		s1,
		wxEmptyString,
		wxEmptyString,
		s3,
		styleFlags
	);
    dialog.CentreOnParent();
//    dialog.SetDirectory(s2);	// screws up getDirectory() in Linux
    if (dialog.ShowModal() == wxID_OK){
        returnDirectory = dialog.GetPath();
        wxFileName filepath( returnDirectory );
        filepath.Normalize(wxPATH_NORM_LONG|wxPATH_NORM_DOTS|wxPATH_NORM_TILDE|wxPATH_NORM_ABSOLUTE);
        filepath.MakeRelativeTo(wxFileName::GetCwd()); 
//        filepath.PrependDir( "." );
        returnDirectory = filepath.GetPathWithSep();
        returnFilename = filepath.GetFullName();
        std::cout << "Return directory: " << returnDirectory << std::endl;
        std::cout << "Return filename: " << returnFilename << std::endl;


/*        filepath.MakeRelativeTo(wxFileName::GetCwd()); 
        filepath.PrependDir( "." );
        returnDirectory = filepath.GetFullPath();

        wxString dir,fnam,fext;
        wxFileName::SplitPath(returnDirectory,&dir,&fnam,&fext);
        std::cout << returnDirectory << std::endl;

        returnDirectory = dir + wxFileName::GetPathSeparator();
        returnFilename = fnam + "." + fext; */
    }
	else	{
		returnDirectory = "";
		returnFilename = "";
	}
}

void MainFrame::setToolState(int toolState)
{
	_toolState = toolState;
	if(toolState==1)
		_toolsMenu->Check(TOOLS_HOOK,true);
	else if(toolState==2)
		_toolsMenu->Check(TOOLS_SCALPEL,true);
	else if(toolState==3)
		_toolsMenu->Check(TOOLS_SUTURE,true);
	else if(toolState==4)
		_toolsMenu->Check(TOOLS_EXCISE,true);
	else
		_toolsMenu->Check(TOOLS_VIEW,true);
	_ftd->setToolBoxSelection(toolState);
	_surgAct->setToolState(toolState);
}

int MainFrame::showModalMessageDialog(const char *message, const char *messageTitle, bool OKnotYesNoButtons, bool cancelButton)
{
	long style = wxICON_INFORMATION|wxSTAY_ON_TOP;
	if(OKnotYesNoButtons)
		style |= wxOK;
	else
		style |= wxYES_NO|wxNO_DEFAULT;
	if(cancelButton)
		style |= wxCANCEL;
	wxString mt1,mt2;
	mt1.append(message);
	mt2.append(messageTitle);
	wxMessageDialog dlg(NULL,mt1,mt2,style);
	int retVal = dlg.ShowModal();
	return retVal;
}

 #if !defined(NO_PHYSICS) && defined(USE_RPC_INTERFACE)
int MainFrame::getRPCServerDialog(std::string& server_address, std::string& server_port)
{
	long style = wxICON_INFORMATION;	// NATHAN - I removed this. |wxSTAY_ON_TOP;
    style |= wxOK;
    style ^= wxTE_PASSWORD;
	
    wxTextEntryDialog dlg(this, "Enter location of remote simulation server: ", "", "http://HOST.cs.wisc.edu:9090", style);
	int retVal = dlg.ShowModal();

    server_address = "localhost";
    server_port = "9090";
    if(retVal == wxID_OK)
        {
            std::cout<<dlg.GetValue()<<std::endl;
            wxURI rpc_server(dlg.GetValue());
            if( rpc_server.HasServer() )
                server_address = rpc_server.GetServer();
            if( rpc_server.HasPort() )
                server_port = rpc_server.GetPort();          
        }
    
    std::cout<<server_address<<std::endl;
    std::cout<<server_port<<std::endl;

	return retVal;
}
#endif

//---------------------------------------------------------------------------
// TestGLCanvas
//---------------------------------------------------------------------------

BEGIN_EVENT_TABLE(BaseGLCanvas, wxGLCanvas)
    EVT_SIZE(BaseGLCanvas::OnSize)
    EVT_PAINT(BaseGLCanvas::OnPaint)
    EVT_KEY_DOWN(BaseGLCanvas::OnKeyDown)
    EVT_KEY_UP(BaseGLCanvas::OnKeyUp)
    EVT_MOUSE_EVENTS(BaseGLCanvas::OnMouseEvent)
    EVT_IDLE(BaseGLCanvas::OnIdle)
END_EVENT_TABLE()

BaseGLCanvas::BaseGLCanvas(wxWindow *parent,
                           wxWindowID id,
                           int* gl_attrib)
    : wxGLCanvas(parent, id, gl_attrib),_glInitialized(false)
{

	_frame = (MainFrame*)parent;
    // Explicitly create a new rendering context instance for this canvas.
    m_glRC = new wxGLContext(this);
	_surgicalDrag = false;
}

BaseGLCanvas::~BaseGLCanvas()
{
    delete m_glRC;
}

void BaseGLCanvas::OnPaint( wxPaintEvent& WXUNUSED(event) )
{
    // This is a dummy, to avoid an endless succession of paint messages.
    // OnPaint handlers must always create a wxPaintDC.
    wxPaintDC dc(this);

    // This is normally only necessary if there is more than one wxGLCanvas
    // or more than one wxGLContext in the application.
//    SetCurrent(*m_glRC);
    if(!_glInitialized)
    {
        SetCurrent(*m_glRC);
        std::string errorMessage;
        if(!_wGG.initializeGraphics(errorMessage))	{
        	_frame->showModalMessageDialog(errorMessage.c_str(),"Base Graphics Error",true,false);
        	_frame->Close();
        }
    	_glInitialized = true;
    }
	_wGG.drawAll();
    SwapBuffers();
}

void BaseGLCanvas::OnSize(wxSizeEvent& event)
{
    // This is normally only necessary if there is more than one wxGLCanvas
    // or more than one wxGLContext in the application.
#ifdef WINDOWS
	SetCurrent(*m_glRC);
#endif
	event.Skip();

    // It's up to the application code to update the OpenGL viewport settings.
    // This is OK here only because there is only one canvas that uses the
    // context. See the cube sample for that case that multiple canvases are
    // made current with one context.
	_xSize = event.GetSize().x;
	_ySize = event.GetSize().y;
	_wGG.setViewport(0,0,_xSize,_ySize);
}

void BaseGLCanvas::OnKeyDown(wxKeyEvent& event)
{
	int key;
    switch( (key = event.GetKeyCode()) )
    {
    case WXK_ESCAPE:
        wxTheApp->ExitMainLoop();
        return;

    case WXK_LEFT:
        break;

    case WXK_RIGHT:
        break;

    case WXK_UP:
        break;

    case WXK_DOWN:
        break;

//    case WXK_DELETE:
//		_surgAct->onKeyDown(WXK_DELETE);
//		break;

    default:
		_surgAct->onKeyDown(key);
    	event.Skip();
    }
    Refresh(false);
    wxIdleEvent iEvent;
    OnIdle(iEvent);
}

void BaseGLCanvas::OnKeyUp(wxKeyEvent& event)
{
    _surgAct->onKeyUp(event.GetKeyCode());
   	event.Skip();
    Refresh(false);
    wxIdleEvent iEvent;
    OnIdle(iEvent);
}

void BaseGLCanvas::OnMouseEvent(wxMouseEvent& event)
{
    // Allow default processing to happen, or else the canvas cannot gain focus
    // (for key events).
	//COURT - remember this!
    event.Skip();

	bool dragging=false;
	int button=-1;
	unsigned short screenX=event.GetX(),screenY=event.GetY();
	
	if(event.Dragging())
		dragging=true;
	if(event.LeftIsDown())
		button=0;
	else if(event.RightDown())	{
		std::string name; float position[3]; int triangle;
		_wGG.pick(screenX,screenY,name,position,triangle);
		if(!name.empty())
		{
			if(_surgAct->rightMouseDown(name,position,triangle))	{
				_lastSurgX = screenX;
				_lastSurgY = screenY;
			       Refresh(false);
				_surgicalDrag = true;
				return;	// surgical action taken.  Don't process this event as a viewer command
			}
		}
	}
	else if(event.RightIsDown())	{
		button=2;
	}
	else if(event.MiddleIsDown())
		button=1;
	else if(event.RightUp())	{
		_surgicalDrag = false;
		std::string name; float position[3]; int triangle;
		_wGG.pick(screenX,screenY,name,position,triangle,true);
		if(_surgAct->rightMouseUp(name,position,triangle))	{
			Refresh(false);
		    wxIdleEvent iEvent;
		    OnIdle(iEvent);
			return;	// Don't process this event as a viewer command
		}
	}
	else
		return;
	if(button==2)	{
		if(dragging) {
			if(_surgicalDrag) {
				_surgAct->mouseMotion((float)(screenX-_lastSurgX)/_xSize,(float)(_lastSurgY-screenY)/_ySize);
				_lastSurgX = screenX;
				_lastSurgY = screenY;
		        Refresh(false);
				return;	// surgical action taken.  Don't process this event as a viewer command 
			}
		}
	}
	_wGG.mouseButtonEvent(screenX,screenY,button,dragging);
	if(dragging)
        Refresh(false);
}

void BaseGLCanvas::OnIdle(wxIdleEvent &event)
{
    if( _surgAct )
        _surgAct->updatePhysics();
    Refresh(false);
	event.RequestMore();
	event.Skip();
}
