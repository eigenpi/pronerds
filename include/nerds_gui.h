#ifndef _GUI_GRAPHICS_H_
#define _GUI_GRAPHICS_H_

#include "nerds.h"
#include "nerds_utils.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>

#define MWIDTH 104 // width of menu window; 
#define T_AREA_HEIGHT 24 // height of text window;
#define MAX_FONT_SIZE 40 // largest point size of text;
#define BUTTON_TEXT_LEN 20
#define MAXPIXEL 15000   
#define MINPIXEL -15000
#define COLORS_NUMBER 11
#define MAXPTS 100 // Maximum number of points drawable by fillpoly;
#define OFF 1
#define ON 0
#define PRIORITY_MINOR 0 // for update_screen; denotes importance of update;
#define PRIORITY_MAJOR 1

#define max(a,b) (((a) > (b))? (a) : (b))
#define min(a,b) ((a) > (b)? (b) : (a))

enum COLOR_TYPES { WHITE, BLACK, DARKGREY, LIGHTGREY, BLUE, GREEN, YELLOW, CYAN,
    RED, GREY35, GREY15 }; // DARKGREEN, MAGENTA
enum LINE_TYPES { SOLID, DASHED };
enum PICTURE_TYPE { NO_PICTURE, NODES }; // What's on screen?

////////////////////////////////////////////////////////////////////////////////
//
// GUI_GRAPHICS
//
////////////////////////////////////////////////////////////////////////////////

class GUI_GRAPHICS {
 public:
    enum DISPLAY_TYPE { SCREEN, POSTSCRIPT };
    enum BUTTON_FUNCTION { UP, DOWN, LEFT, RIGHT, ZOOM_IN, ZOOM_OUT, ZOOM_FIT,
        ADJUSTWIN, DRAWSCREEN, PS, PROCEED, TOGGLE_TIES, 
        QUIT };
 public:
    class BUTTON_TYPE {
    public:
        int width; 
        int height; 
        int xleft; 
        int ytop;
        BUTTON_FUNCTION fcn;
        Window win; 
        int istext; 
        char text[BUTTON_TEXT_LEN]; 
        int ispoly; 
        int poly[3][2]; 
        int ispressed;
    public:
        BUTTON_TYPE() {}
        ~BUTTON_TYPE() {}           
    };

    typedef struct { float x; float y; } T_POINT; // Used in calls to fillpoly

 private:
    bool _is_gui_usable;
    int _menu_font_size; // font for menus and dialog boxes;
    int _num_buttons; // Number of menu buttons
    BUTTON_TYPE *_button; // [0..num_buttons-1]
    DISPLAY_TYPE _display_type; // for SCREEN or POSTSCRIPT?
    int _screen_num;
    GC _gc;
    GC _gcxor; 
    GC _gc_menus;
    XFontStruct *_font_info[MAX_FONT_SIZE + 1]; // Data for each size 
    int _font_is_loaded[MAX_FONT_SIZE + 1]; // 1: loaded, 0: not
    unsigned int _display_width;
    unsigned int _display_height; // screen size;
    unsigned int _top_width;
    unsigned int _top_height; // window size;
    Window _toplevel;
    Window _menu;
    Window _textarea; // various windows;
    // world coordinates
    float _xleft, _xright, _ytop, _ybot; 
    // Initial world coordinates
    float _saved_xleft, _saved_xright, _saved_ytop, _saved_ybot; 
    // Figure boundaries for PostScript output, in PostScript coordinates
    float _ps_left, _ps_right, _ps_top, _ps_bot;
    float _ps_xmult, _ps_ymult; // Transformation for PostScript.
    float _xmult, _ymult; // Transformation factors
    Colormap _private_cmap; // "None" unless a private cmap was allocated;
    // Graphics state.  Set start-up defaults here.
    int _currentcolor;
    int _currentlinestyle;
    int _currentlinewidth;
    int _currentfontsize;
    char _user_message[BUFFER_SIZE]; // User message to display
    int _colors[COLORS_NUMBER]; // Color indices passed back from X Windows;
    // For PostScript output
    FILE *_ps_file;
    int _piccount; // keeps track of saved .ps pictures;

    char _default_message[BUFFER_SIZE]; // Default screen message on screen
    PICTURE_TYPE _pic_on_screen;
    bool _show_ties; // display tie switches?

    // sketch storage arrays for speed up of gui;
    // the whole system of coordinates is for now like a matrix of tiles
    // where nodes may be placed; this is to have a regularity; these two
    // arrays store pre-computed left,bottom coordinates to minimize the
    // amount of on-the-fly calculations;
    // COORDINATE SYSTEM goes from (0,0) at the lower left corner to
    // (x_tile_left[nx+1] + tile_width, y_tile_bottom[ny+1] + tile_width) 
    // in the upper right corner.
    vector<float> _x_tile_left;
    vector<float> _y_tile_bottom;
    vector<int> _x_coords;
    vector<int> _y_coords;
    vector<COLOR_TYPES> _node_color;
    vector<double> _node_voltage;
    vector<int> _x_src;
    vector<int> _y_src;
    vector<int> _x_des;
    vector<int> _y_des;
    vector<COLOR_TYPES> _arc_color; 
    vector<LINE_TYPES> _arc_line;
    // numerical values controlled from main function; setting the nodes size,
    // used during drawing;
    float _tile_width;
    float _channel_width;
    // Need user input after (such as proceed): 
    // 0: each time unit, 1: each print-results cycle, 2: never
    int _wait_for_user_input_automode;
 
 public:
    Display *_display; // need it public as it's given to X routines;

 public:
    PRONERDS *_pronerds; // the main program platform;


 public:
    GUI_GRAPHICS( PRONERDS *pronerds);
    ~GUI_GRAPHICS() {}


    void set_is_gui_usable(bool val) { _is_gui_usable = val; }
    bool is_gui_usable() { return _is_gui_usable; }
    bool build();
    void build_textarea(void);
    void build_default_menu (void);
    void setpoly(int bnum, int xc, int yc, int r, float theta);
    void turn_on_off (int pressed);
    int which_button (Window win);  
    void drawmenu(void);
    void update_transform (void);
    void update_ps_transform (void);
    int rect_off_screen (float x1, float y1, float x2, float y2);

    void set_graphics_state( bool show_graphics_val, int wfui_automode_val);
    void update_screen( int priority, char *msg,
        PICTURE_TYPE pic_on_screen_val);
    void drawscreen (void);
    void redraw_screen (void);
    void toggle_ties(void);
    void alloc_draw_structs (void);
    void deselect_all (void);
    void init_draw_coords (float tile_width_val, float channel_width_val);
    void draw_nodes(void);
    void draw_ties(void);

    // Routine for X Windows Input.  act_on_button responds to buttons 
    // being pressed in the graphics area.  drawscreen is the user's   
    // routine that can redraw all the graphics.
    void event_loop( void);
    void execute_action_associated_with_button_function( BUTTON_FUNCTION fcn);
    void create_button (char *prev_button_text, char *button_text, 
        BUTTON_FUNCTION button_func);
    void destroy_button (char *button_text);


    void init_graphics (char *window_name); // Initializes X display  
    void close_graphics (void); // Closes X display
    // Changes message in text area. 
    void update_message (char *msg);
    // Normal users shouldn't have to use draw_message.  Should only be
    // useful if using non-interactive graphics and you want to redraw 
    // yourself because of an expose.
    void draw_message (void);
    // Sets world coordinates 
    void init_world (float xl, float yt, float xr, float yb);
    void flushinput (void); // Empties event queue

    int xcoord(float worldx);
    int ycoord(float worldy);
    void load_font(int pointsize);
    void force_setcolor(int cindex);
    void force_setlinestyle(int linestyle);
    void force_setlinewidth (int linewidth);
    void force_setfontsize (int pointsize);
    // Following routines draw to SCREEN if disp_type = SCREEN
    // and to a PostScript file if disp_type = POSTSCRIPT
    void setcolor (int cindex); // Use a constant from clist 
    void setlinestyle (int linestyle);  
    void setlinewidth (int linewidth);
    void setfontsize (int pointsize);
    void drawline (float x1, float y1, float x2, float y2);
    void drawrect (float x1, float y1, float x2, float y2);
    void fillrect (float x1, float y1, float x2, float y2);
    void fillpoly (T_POINT *points, int npoints);
    // Draw or fill a circular arc, respectively.  Angles in degrees.  startang 
    // measured from positive x-axis of Window.  Positive angextent means       
    // counterclockwise arc.
    void drawarc (float xcen, float ycen, float rad, float startang, float angextent); 
    void fillarc (float xcen, float ycen, float rad, float startang, float angextent);

    float angnorm (float ang);
    // boundx specifies horizontal bounding box.  If text won't fit in   
    // the space specified by boundx (world coordinates) text isn't drawn
    void drawtext (float xc, float yc, char *text, float boundx);
    void clearscreen (void); // Erases the screen 

    // Opens file for postscript commands and initializes it.  All subsequent 
    // drawing commands go to this file until close_postscript is called.
    int init_postscript (char *fname); // Returns 1 if successful
    // Closes file and directs output to screen again.
    void close_postscript (void);

    // Function declarations for button responses. 
    void update_win (int x[2], int y[2]);
    void translate_up (); 
    void translate_left (); 
    void translate_right (); 
    void translate_down (); 
    void zoom_in ();
    void zoom_out ();
    void zoom_fit ();
    void adjustwin (); 
    void postscript ();
    void proceed ();
    void quit ();

    //Bool test_if_exposed (Display *disp, XEvent *event_ptr, XPointer dummy);
    inline void map_button (int bnum);
    inline void unmap_button (int bnum);

    void menutext(Window win, int xc, int yc, char *text);
    void drawbut (int bnum);

    // Macros for translation from world to PostScript coordinates
    inline float XPOST( float worldx)
        { return (((worldx) - _xleft) * _ps_xmult + _ps_left); }
    inline float YPOST( float worldy) 
        { return (((worldy) - _ybot) * _ps_ymult + _ps_bot); }
    // Macros to convert from X Windows Internal Coordinates to my
    // World Coordinates.  (This macro is used only rarely, so
    // the divides don't hurt speed).
    inline float XTOWORLD(int x) { return (((float) x) / _xmult + _xleft); }
    inline float YTOWORLD(int y) { return (((float) y) / _ymult + _ytop); }
};

#endif
