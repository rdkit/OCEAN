// dragable object: see http://aktuell.de.selfhtml.org/artikel/javascript/draganddrop/index.htm
var dragobjekt = null;

var DBI_dragx = 0;
//var DBI_dragy = 0;
var DBI_posx = 0;
//var DBI_posy = 0;

var DBI_sizeX = 0;
//var DBI_sizeY = 0;

function draginit() {
 // Initialisierung der überwachung der Events
  document.onmousemove = drag;
  document.onmouseup = dragstop;
}
function dragstart(element) {
   //Wird aufgerufen, wenn ein Objekt bewegt werden soll.

  dragobjekt = element;
  DBI_dragx = DBI_posx - dragobjekt.offsetLeft;
//  DBI_dragy = DBI_posy - dragobjekt.offsetTop;
  DBI_sizeX = dragobjekt.width;
//  DBI_sizeY = dragobjekt.height;
}

function dragstop() {
  //Wird aufgerufen, wenn ein Objekt nicht mehr bewegt werden soll.

//  alert(dragobjekt);
  dragobjekt=null;

}

function drag(ereignis) {
  //Wird aufgerufen, wenn die Maus bewegt wird und bewegt bei Bedarf das Objekt.

  DBI_posx = document.all ? window.event.clientX : ereignis.pageX;
//  DBI_posy = document.all ? window.event.clientY : ereignis.pageY;
  if(dragobjekt != null) {
    dragobjekt.style.left = (DBI_posx - DBI_dragx) + "px";
//    dragobjekt.style.top = (DBI_posy - DBI_dragy) + "px";
//    if (dragobjekt.style.top.startsWith("-")){
//        dragobjekt.style.top = 0 + "px";
//        console.log("ds " + dragobjekt.style.top);
//    }
//    console.log("ds " + dragobjekt.style.top);
//    dragobjekt.width = DBI_sizeX;
//    dragobjekt.height = DBI_sizeY;
	new_width = $(dragobjekt).width();
	new_height = $(dragobjekt).height();
	elems = dragobjekt.getElementsByTagName("img");
    for (i=0;i<elems.length;i++){
        percentage = new_width / DBI_sizeX;
        $(elems[i]).width(new_width - 20);
//        alert(new_width);
        plot_pie(new_width,new_width)
//        console.log((new_width) + " ; " + (new_width));
//	    console.log("img.width "+elems[i].width);
//	    console.log("new.width "+new_width);
//	    console.log("new.height "+new_height);
//	    console.log("percentage " + percentage);
    }
  }

}

