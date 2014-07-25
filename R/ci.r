
exactci = function(predictions, observations){
	out = .Call("exactciR2C", predictions, observations)
	return (out); 
}

equalci = function(predictions, observations){
	out = .Call("equalciR2C", predictions, observations)
	return (out); 
}
