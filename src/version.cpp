//
// Created by heumos on 01.12.21.
//

// modified from https://github.com/vgteam/vg/blob/master/src/version.cpp

#include "version.hpp"

// Get the git version macro from the build system
#include "../include/smoothxg_git_version.hpp"

#include <iostream>
#include <sstream>

// If the smoothxg_GIT_VERSION deosn't exist at all, define a placeholder
// This lets us be somewhat robust to undeterminable versions
#ifndef SMOOTHXG_GIT_VERSION
#define SMOOTHXG_GIT_VERSION "not-from-git"
#endif

// Define a way to quote macro values.
// See https://stackoverflow.com/a/196093
//#define QUOTE(arg) #arg
// We need another level to get the macro's value and not its name.
//#define STR(macro) QUOTE(macro)

namespace smoothxg {

	using namespace std;

	// Define all the strings as the macros' values
	const string Version::VERSION = SMOOTHXG_GIT_VERSION;

	// Keep the list of codenames.
	// Add new codenames here
	const unordered_map<string, string> Version::codenames = {
			{"v0.1", "initial release"},
			{"v0.2", "edit-distance based block splitting"},
			{"v0.3-alpha", "super smoothxg pre-release"},
			{"v0.4", "smoothxg 0.4 - super"},
			{"v0.5", "stabilizing the smooth"},
			{"v0.5.1", "Smoothxg 0.5.1 - Scaling up"},
			{"v0.5.2", "Smoothxg 0.5.2 - Imbottitura"},
			{"v0.6.0", "unciampato"},
			{"v0.6.1", "Smoothxg 0.6.1 - Magro"},
			{"v0.6.2", "Magrissimo"},
			{"v0.6.3", "Generico"},
			{"v0.6.4", "Pasticcione"}
			// Add more codenames here
	};

	string Version::get_version() {
		return VERSION;
	}

	string Version::get_release() {
		auto dash = VERSION.find('-');
		if (dash == -1) {
			// Pure tag versions have no dash
			return VERSION;
		} else {
			// Otherwise it is tag-count-hash and the tag describes the release
			return VERSION.substr(0, dash);
		}
	}

	string Version::get_codename() {
		auto release = get_release();

		auto found = codenames.find(release);

		if (found == codenames.end()) {
			// No known codename for this release.
			// Return an empty string so we can just not show it.
			return "";
		} else {
			// We have a known codename!
			return found->second;
		}
	}

	string Version::get_short() {
		stringstream s;
		s << VERSION;

		auto codename = get_codename();
		if (!codename.empty()) {
			// Add the codename if we have one
			s << " \"" << codename << "\"";
		}

		return s.str();
	}

}
