#include "color.h"

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
            fprintf(stderr, "ℹ [INFO] ");
            break;
        case DEBUG_SUCCESS:
            ncdirect_set_fg_rgb(ctx->nd, THEME_SUCCESS);
            fprintf(stderr, "✓ [OK] ");
            break;
        case DEBUG_WARNING:
            ncdirect_set_fg_rgb(ctx->nd, THEME_WARNING);
            fprintf(stderr, "⚠ [WARN] ");
            break;
        case DEBUG_ERROR:
            ncdirect_set_fg_rgb(ctx->nd, THEME_ERROR);
            fprintf(stderr, "✗ [ERROR] ");
            break;
        case DEBUG_TRACE:
            ncdirect_set_fg_rgb(ctx->nd, THEME_TRACE);
            fprintf(stderr, "→ [TRACE] ");
            break;
    }
    
    ncdirect_set_fg_rgb(ctx->nd, THEME_TEXT);
    vfprintf(stderr, fmt, args);
    fprintf(stderr, "\n");
    
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
    fprintf(stderr, "\n╔");
    for (int i = 0; i < ctx->header_width - 2; i++) fprintf(stderr, "═");
    fprintf(stderr, "╗\n");
    
    // Title with padding
    fprintf(stderr, "║ %-*s ║\n", ctx->header_width - 4, title);
    
    // Bottom border
    fprintf(stderr, "╚");
    for (int i = 0; i < ctx->header_width - 2; i++) fprintf(stderr, "═");
    fprintf(stderr, "╝\n\n");
}

void debug_variable(debug_ctx_t* ctx, const char* name, const char* value) {
    ncdirect_set_fg_rgb(ctx->nd, THEME_VAR_NAME);
    fprintf(stderr, "  %s: ", name);
    ncdirect_set_fg_rgb(ctx->nd, THEME_VAR_VALUE);
    fprintf(stderr, "%s\n", value);
}

void debug_cleanup(debug_ctx_t* ctx) {
    ncdirect_stop(ctx->nd);
    free(ctx);
}

// Example usage


// Compile with: gcc -o debug debug.c -lnotcurses-core
