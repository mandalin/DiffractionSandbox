//
//  CookieExample.js
//  Tangle
//
//  Created by Bret Victor on 6/10/11.
//  (c) 2011 Bret Victor.  MIT open-source license.
//

window.addEvent('domready', function () {

    var model = {
        initialize: function () {
            this.cookies = 3;
			this.c = 343;
			this.r_R=1;
            this.r_S=1;
			this.z_R=1;
            this.z_S=1;
			this.theta_R=1;
            this.theta_S=1;
			this.eta=1;
            this.tau=10;
        },
        update: function () {
            //this.calories = this.cookies * this.caloriesPerCookie;
            //this.dailyPercent = 100 * this.calories / this.caloriesPerDay;
        }
    };
    
    for (var i = 1; ; i++) {
        var id = "SvenssonExample" + ((i > 1) ? i : ""); //ternary
        var element = document.getElementById(id);
        if (!element) { break; }
        new Tangle(element,model);
    }
    
});
