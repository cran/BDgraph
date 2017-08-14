## ----------------------------------------------------------------------------|
# For discrete data: to transfer raw data for the algorithm 
## ----------------------------------------------------------------------------|
transfer = function( r_data )
{
	if( class( r_data ) == "sim" ) r_data <- r_data $ data
  
#	all_patterns    = apply( r_data, 1, function( x ){ paste( x, collapse = '' ) } )   
#	unique_patterns = unique( all_patterns )
	   
#	length_unique_patterns = length( unique_patterns )
#	data  = matrix( 0, nrow = length_unique_patterns, ncol = ncol( r_data ) + 1 )
	   
#	for( i in seq_len( length_unique_patterns ) )
#	{
#		which_one = which( all_patterns == unique_patterns[i] )
#		data[i, ] = c( r_data[ which_one[1], ], length( which_one ) )
#	}

	n    = dim( r_data )[1]
	p    = dim( r_data )[2]
	data = matrix( 0, nrow = n, ncol = p + 1 )
	size_unique_data = 0
	
	result = .C( "transfer_data", as.integer(r_data), data = as.integer(data), as.integer(n), as.integer(p), size_unique_data = as.integer(size_unique_data), PACKAGE = "BDgraph" )
	
	size_unique_data = result $ size_unique_data
	
	label = colnames( r_data )
	if( is.null( label ) ) label = 1:p
  
	data = matrix ( result $ data, n, p + 1, dimnames = list( NULL, c( label, "ferq" ) ) ) 
	data = data[ (1:size_unique_data), ]

	return( data )
}
   
