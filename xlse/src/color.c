#include <notcurses/notcurses.h>
#include <notcurses/direct.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

typedef enum {
    DEBUG_INFO,
    DEBUG_SUCCESS,
    DEBUG_WARNING,
    DEBUG_ERROR,
    DEBUG_TRACE
} debug_level_t;

// Color Themes
// Choose one by uncommenting the desired theme

// ===== THEME: DEFAULT (Vibrant) =====
#define THEME_INFO      0x00bfff  // Sky blue
#define THEME_SUCCESS   0x00ff00  // Green
#define THEME_WARNING   0xffaa00  // Orange
#define THEME_ERROR     0xff0000  // Red
#define THEME_TRACE     0x888888  // Gray
#define THEME_HEADER    0xff00ff  // Magenta
#define THEME_VAR_NAME  0x00ffff  // Cyan
#define THEME_VAR_VALUE 0xffff00  // Yellow
#define THEME_TEXT      0xffffff  // White

/* ===== THEME: GRUVBOX =====
#define THEME_INFO      0x83a598  // Blue
#define THEME_SUCCESS   0xb8bb26  // Green
#define THEME_WARNING   0xfabd2f  // Yellow
#define THEME_ERROR     0xfb4934  // Red
#define THEME_TRACE     0x928374  // Gray
#define THEME_HEADER    0xd3869b  // Purple
#define THEME_VAR_NAME  0x8ec07c  // Aqua
#define THEME_VAR_VALUE 0xfe8019  // Orange
#define THEME_TEXT      0xebdbb2  // Fg
*/

/* ===== THEME: NORD =====
#define THEME_INFO      0x81a1c1  // Blue
#define THEME_SUCCESS   0xa3be8c  // Green
#define THEME_WARNING   0xebcb8b  // Yellow
#define THEME_ERROR     0xbf616a  // Red
#define THEME_TRACE     0x4c566a  // Gray
#define THEME_HEADER    0xb48ead  // Purple
#define THEME_VAR_NAME  0x88c0d0  // Frost
#define THEME_VAR_VALUE 0xd08770  // Orange
#define THEME_TEXT      0xeceff4  // Snow
*/

/* ===== THEME: DRACULA =====
#define THEME_INFO      0x8be9fd  // Cyan
#define THEME_SUCCESS   0x50fa7b  // Green
#define THEME_WARNING   0xf1fa8c  // Yellow
#define THEME_ERROR     0xff5555  // Red
#define THEME_TRACE     0x6272a4  // Comment
#define THEME_HEADER    0xff79c6  // Pink
#define THEME_VAR_NAME  0xbd93f9  // Purple
#define THEME_VAR_VALUE 0xffb86c  // Orange
#define THEME_TEXT      0xf8f8f2  // Foreground
*/

/* ===== THEME: SOLARIZED DARK =====
#define THEME_INFO      0x268bd2  // Blue
#define THEME_SUCCESS   0x859900  // Green
#define THEME_WARNING   0xb58900  // Yellow
#define THEME_ERROR     0xdc322f  // Red
#define THEME_TRACE     0x586e75  // Base01
#define THEME_HEADER    0xd33682  // Magenta
#define THEME_VAR_NAME  0x2aa198  // Cyan
#define THEME_VAR_VALUE 0xcb4b16  // Orange
#define THEME_TEXT      0x839496  // Base0
*/

/* ===== THEME: TOKYO NIGHT =====
#define THEME_INFO      0x7aa2f7  // Blue
#define THEME_SUCCESS   0x9ece6a  // Green
#define THEME_WARNING   0xe0af68  // Yellow
#define THEME_ERROR     0xf7768e  // Red
#define THEME_TRACE     0x565f89  // Comment
#define THEME_HEADER    0xbb9af7  // Purple
#define THEME_VAR_NAME  0x7dcfff  // Cyan
#define THEME_VAR_VALUE 0xff9e64  // Orange
#define THEME_TEXT      0xc0caf5  // Foreground
*/

/* ===== THEME: CATPPUCCIN MOCHA =====
#define THEME_INFO      0x89b4fa  // Blue
#define THEME_SUCCESS   0xa6e3a1  // Green
#define THEME_WARNING   0xf9e2af  // Yellow
#define THEME_ERROR     0xf38ba8  // Red
#define THEME_TRACE     0x6c7086  // Overlay0
#define THEME_HEADER    0xcba6f7  // Mauve
#define THEME_VAR_NAME  0x94e2d5  // Teal
#define THEME_VAR_VALUE 0xfab387  // Peach
#define THEME_TEXT      0xcdd6f4  // Text
*/

typedef struct {
    struct ncdirect* nd;
    int header_width;
} debug_ctx_t;

debug_ctx_t* debug_init(void) {
    debug_ctx_t* ctx = malloc(sizeof(debug_ctx_t));
    ctx->nd = ncdirect_init(NULL, NULL, 0);
    ctx->header_width = -1; // Not yet determined
    return ctx;
}

void debug_log(debug_ctx_t* ctx, debug_level_t level, const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
    
    switch(level) {
        case DEBUG_INFO:
            ncdirect_set_fg_rgb(ctx->nd, THEME_INFO);
            printf("ℹ [INFO] ");
            break;
        case DEBUG_SUCCESS:
            ncdirect_set_fg_rgb(ctx->nd, THEME_SUCCESS);
            printf("✓ [OK] ");
            break;
        case DEBUG_WARNING:
            ncdirect_set_fg_rgb(ctx->nd, THEME_WARNING);
            printf("⚠ [WARN] ");
            break;
        case DEBUG_ERROR:
            ncdirect_set_fg_rgb(ctx->nd, THEME_ERROR);
            printf("✗ [ERROR] ");
            break;
        case DEBUG_TRACE:
            ncdirect_set_fg_rgb(ctx->nd, THEME_TRACE);
            printf("→ [TRACE] ");
            break;
    }
    
    ncdirect_set_fg_rgb(ctx->nd, THEME_TEXT);
    vprintf(fmt, args);
    printf("\n");
    
    va_end(args);
}

void debug_header(debug_ctx_t* ctx, const char* title) {
    int title_len = strlen(title);
    int term_width = ncdirect_dim_x(ctx->nd);
    
    // Determine header width on first call
    if (ctx->header_width == -1) {
        // Use minimum of: terminal width - 4, or title length + 10, capped at 80
        int preferred = title_len + 10;
        int max_width = term_width - 4;
        ctx->header_width = (preferred < max_width) ? preferred : max_width;
        if (ctx->header_width > 80) ctx->header_width = 80;
        if (ctx->header_width < 20) ctx->header_width = 20;
    }
    
    // Ensure we have enough space for the title
    int content_width = ctx->header_width - 4; // Account for "║ " and " ║"
    if (title_len > content_width) {
        // Truncate title if too long
        char truncated[content_width + 1];
        snprintf(truncated, content_width - 2, "%s", title);
        strcpy(truncated + content_width - 3, "...");
        title = truncated;
        title_len = content_width;
    }
    
    ncdirect_set_fg_rgb(ctx->nd, THEME_HEADER);
    
    // Top border
    printf("\n╔");
    for (int i = 0; i < ctx->header_width - 2; i++) printf("═");
    printf("╗\n");
    
    // Title with padding
    printf("║ %-*s ║\n", ctx->header_width - 4, title);
    
    // Bottom border
    printf("╚");
    for (int i = 0; i < ctx->header_width - 2; i++) printf("═");
    printf("╝\n\n");
}

void debug_variable(debug_ctx_t* ctx, const char* name, const char* value) {
    ncdirect_set_fg_rgb(ctx->nd, THEME_VAR_NAME);
    printf("  %s: ", name);
    ncdirect_set_fg_rgb(ctx->nd, THEME_VAR_VALUE);
    printf("%s\n", value);
}

void debug_cleanup(debug_ctx_t* ctx) {
    ncdirect_stop(ctx->nd);
    free(ctx);
}

// Example usage
int main(void) {
    debug_ctx_t* debug = debug_init();
    
    debug_header(debug, "Application Startup");
    
    debug_log(debug, DEBUG_INFO, "Initializing system...");
    debug_log(debug, DEBUG_TRACE, "Loading configuration from config.json");
    debug_log(debug, DEBUG_SUCCESS, "Configuration loaded successfully");
    
    debug_header(debug, "Processing Data");
    
    debug_variable(debug, "Status", "Active");
    debug_variable(debug, "Records", "1,234");
    debug_variable(debug, "Memory", "45.2 MB");
    
    debug_log(debug, DEBUG_WARNING, "Cache is 80% full");
    debug_log(debug, DEBUG_ERROR, "Failed to connect to database");
    debug_log(debug, DEBUG_INFO, "Retrying connection...");
    debug_log(debug, DEBUG_SUCCESS, "Connected to database");
    
    debug_header(debug, "Shutdown");
    debug_log(debug, DEBUG_INFO, "Cleaning up resources...");
    
    debug_cleanup(debug);
    return 0;
}

// Compile with: gcc -o debug debug.c -lnotcurses-core
