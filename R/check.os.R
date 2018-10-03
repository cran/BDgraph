check.os <- function( os = 2 )
{
	tmp <- .C( 'check_os', os = as.integer( os ), PACKAGE = "BDgraph" )
	
	if( !( tmp $ os %in% 0:1 ) ) stop( 'Failed to indentify the OS. Please contact the authors. ' )
	
	systm <- ifelse( tmp $ os == 0, 'windows_or_mac', 'linux' )
	message( paste0( '   This package ', ifelse( systm == 'windows_or_mac', 'does not support ', 'supports ' ), 'multi-threading on this OS' ) )
}
   
