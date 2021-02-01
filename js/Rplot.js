(function($) {
    $(document).ready(function() {
	
	$('#Rplot').scianimator({
	    'images': ['images/Rplot1.png', 'images/Rplot2.png'],
	    'width': 480,
	    'delay': 1000,
	    'loopMode': 'none',
  'controls': ['first', 'previous', 
         'play', 'next', 'last', 'loop',
         'speed'], 'delayMin' : 0
	});
    });
})(jQuery);
