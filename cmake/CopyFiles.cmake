macro(CopyDLL target_name)
    if (WIN32)
        add_custom_command(
            TARGET ${target_name} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy ${SDL2_ROOT}/SDL2.dll $<TARGET_FILE_DIR:${target_name}>)
        # If you have selected SDL2 component when installed Vulkan SDK, the command as follows will work
        # add_custom_command(
        #     TARGET ${target_name} POST_BUILD
        #     COMMAND ${CMAKE_COMMAND} -E copy ${SDL2_BIN_DIR}/SDL2.dll $<TARGET_FILE_DIR:${target_name}>)
    endif()
endmacro(CopyDLL)