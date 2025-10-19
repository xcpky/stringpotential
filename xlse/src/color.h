#ifndef COLOR_H
#define COLOR_H

#include <stdarg.h>

#define GRUVBOX

typedef enum {
    DEBUG_INFO,
    DEBUG_SUCCESS,
    DEBUG_WARNING,
    DEBUG_ERROR,
    DEBUG_TRACE
} debug_level_t;

// Color Themes
// Choose one by uncommenting the desired theme

#ifdef VIBRANT
#define THEME_INFO 0x00bfff      // Sky blue
#define THEME_SUCCESS 0x00ff00   // Green
#define THEME_WARNING 0xffaa00   // Orange
#define THEME_ERROR 0xff0000     // Red
#define THEME_TRACE 0x888888     // Gray
#define THEME_HEADER 0xff00ff    // Magenta
#define THEME_VAR_NAME 0x00ffff  // Cyan
#define THEME_VAR_VALUE 0xffff00 // Yellow
#define THEME_TEXT 0xffffff      // White
#elif defined(GRUVBOX)

#define THEME_INFO 0x83a598      // Blue
#define THEME_SUCCESS 0xb8bb26   // Green
#define THEME_WARNING 0xfabd2f   // Yellow
#define THEME_ERROR 0xfb4934     // Red
#define THEME_TRACE 0x928374     // Gray
#define THEME_HEADER 0xd3869b    // Purple
#define THEME_VAR_NAME 0x8ec07c  // Aqua
#define THEME_VAR_VALUE 0xfe8019 // Orange
#define THEME_TEXT 0xebdbb2      // Fg

#elif defined(NORD)

#define THEME_INFO 0x81a1c1      // Blue
#define THEME_SUCCESS 0xa3be8c   // Green
#define THEME_WARNING 0xebcb8b   // Yellow
#define THEME_ERROR 0xbf616a     // Red
#define THEME_TRACE 0x4c566a     // Gray
#define THEME_HEADER 0xb48ead    // Purple
#define THEME_VAR_NAME 0x88c0d0  // Frost
#define THEME_VAR_VALUE 0xd08770 // Orange
#define THEME_TEXT 0xeceff4      // Snow

#elif defined(DRACULA)

#define THEME_INFO 0x8be9fd      // Cyan
#define THEME_SUCCESS 0x50fa7b   // Green
#define THEME_WARNING 0xf1fa8c   // Yellow
#define THEME_ERROR 0xff5555     // Red
#define THEME_TRACE 0x6272a4     // Comment
#define THEME_HEADER 0xff79c6    // Pink
#define THEME_VAR_NAME 0xbd93f9  // Purple
#define THEME_VAR_VALUE 0xffb86c // Orange
#define THEME_TEXT 0xf8f8f2      // Foreground

#elif defined(SOLARIZED_DARK)
#define THEME_INFO 0x268bd2      // Blue
#define THEME_SUCCESS 0x859900   // Green
#define THEME_WARNING 0xb58900   // Yellow
#define THEME_ERROR 0xdc322f     // Red
#define THEME_TRACE 0x586e75     // Base01
#define THEME_HEADER 0xd33682    // Magenta
#define THEME_VAR_NAME 0x2aa198  // Cyan
#define THEME_VAR_VALUE 0xcb4b16 // Orange
#define THEME_TEXT 0x839496      // Base0

#elif defined(TOKYO_NIGHT)
#define THEME_INFO 0x7aa2f7      // Blue
#define THEME_SUCCESS 0x9ece6a   // Green
#define THEME_WARNING 0xe0af68   // Yellow
#define THEME_ERROR 0xf7768e     // Red
#define THEME_TRACE 0x565f89     // Comment
#define THEME_HEADER 0xbb9af7    // Purple
#define THEME_VAR_NAME 0x7dcfff  // Cyan
#define THEME_VAR_VALUE 0xff9e64 // Orange
#define THEME_TEXT 0xc0caf5      // Foreground
#elif defined(CATPPUCCIN_MOCHA)

#define THEME_INFO 0x89b4fa      // Blue
#define THEME_SUCCESS 0xa6e3a1   // Green
#define THEME_WARNING 0xf9e2af   // Yellow
#define THEME_ERROR 0xf38ba8     // Red
#define THEME_TRACE 0x6c7086     // Overlay0
#define THEME_HEADER 0xcba6f7    // Mauve
#define THEME_VAR_NAME 0x94e2d5  // Teal
#define THEME_VAR_VALUE 0xfab387 // Peach
#define THEME_TEXT 0xcdd6f4      // Text
#endif

typedef struct {
    struct ncdirect *nd;
    int header_width;
} debug_ctx_t;

debug_ctx_t *debug_init(void);
void debug_log(debug_ctx_t *, debug_level_t, const char *fmt, ...);
void debug_header(debug_ctx_t*, const char*);
void debug_variable(debug_ctx_t*, const char* name, const char*value);
void debug_cleanup(debug_ctx_t*);

#endif
